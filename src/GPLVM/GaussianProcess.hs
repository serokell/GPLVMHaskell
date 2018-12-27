{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE ConstraintKinds #-}

module GPLVM.GaussianProcess
       ( GaussianProcess (..)
       , PosteriorSample (..)
       , functionalPrior
       , gpToTrainingData
       , kernelGP
       , meanGP
       , testValueCovariaceMatrix
       ) where

import Universum hiding (transpose)

import Control.Lens (makeLenses)

import Data.Array.Repa
import Data.Array.Repa.Repr.Unboxed (Unbox)
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import System.Random (Random, RandomGen)

import GPLVM.Types (KernelFunction, Matrix (..), Mean,
                    InputObservations (..), unInputObs, unMatrix)
import GPLVM.Util (randomMatrixD)

data GaussianProcess a = GaussianProcess
    { _kernelGP :: KernelFunction a
    , _meanGP ::  Mean D a
    }

data GPTrainingData a = GPTrainingData
    { _inputTrain :: Matrix D a
    , _outputTrain :: Matrix D a
    }

makeLenses ''GaussianProcess
makeLenses ''GPTrainingData

newtype PosteriorSample a = PosteriorSample (Matrix D a)

-- covariance between gP input and gP output
-- (just apply kernel function)

testValueCovariaceMatrix
    :: forall a. GaussianProcess a
    -> InputObservations a
    -> Matrix D a
testValueCovariaceMatrix gP observations = kernelMatrix
    where
        kernelMatrix = kernelFun points points
        kernelFun = gP ^. kernelGP
        points = observations ^. unInputObs

type GPConstraint a =
    ( Field a
    , Random a
    , Unbox a
    , Floating a
    )

functionalPrior
    :: forall a g.
    ( GPConstraint a
    , RandomGen g
    )
    => Matrix D a
    -> g
    -> Maybe (Matrix D a)
functionalPrior matrix@(Matrix (ADelayed (Z :. rows :. cols) _)) gen = do
    case (mbChol . symS) $ matrix ^. unMatrix of
        Nothing -> Nothing
        Just matrix' -> do
            let matrixD = smap id matrix'
            return $ Matrix $ smap id (matrixD `mulS` (randomCoeffs ^. unMatrix))
    where
        randomCoeffs = randomMatrixD gen (cols, rows)

gpToTrainingData
    :: forall a g.
    ( GPConstraint a
    , Eq a
    , RandomGen g
    )
    => InputObservations a
    -> GaussianProcess a
    -> GPTrainingData a
    -> g
    -> Maybe (PosteriorSample a)
gpToTrainingData inputObserve gP trainingData gen = do
    let covarianceMatrix = testValueCovariaceMatrix gP inputObserve
        -- get covariance matrix for test inputs
    let inputTrain' = trainingData ^. inputTrain

        -- kernel applied to input training points
    let trainingKernel = (gP ^. kernelGP) inputTrain' inputTrain'
        -- Cholesky decomposition applied to kernel of training points
    let cholK = (cholD . symS) (trainingKernel ^. unMatrix)

        -- test points mean
    let testPointMean = (gP ^. kernelGP) (inputObserve ^. unInputObs) inputTrain'

    cholKSolve <- linearSolveS cholK (testPointMean ^. unMatrix)

        -- solve linear system for output training points
    cholKSolveOut <- linearSolveS cholK ((trainingData ^. outputTrain) ^. unMatrix)

    let cholKSolve' = smap id cholKSolve
    let cholKSolveOut' = smap id cholKSolveOut
    let mulD m n = smap id (mulS m n)

        -- compute mean
    let mean = transpose cholKSolve' `mulD` cholKSolveOut'

        -- posterior
    let postF' = (cholD . symS) $
                     (covarianceMatrix ^. unMatrix) -^
                     (transpose cholKSolve' `mulD` cholKSolve')
    prior' <- functionalPrior (Matrix postF') gen
    let prior = prior' ^. unMatrix
    return $ PosteriorSample . Matrix $ mean +^ prior
