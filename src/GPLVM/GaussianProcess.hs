{-# LANGUAGE AllowAmbiguousTypes, ConstraintKinds #-}

module GPLVM.GaussianProcess
       ( GaussianProcess (..)
       , GPMean
       , GPTrainingData (..)
       , PosteriorSample (..)
       , functionalPrior
       , gpToPosteriorSample
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

import GPLVM.Types (InputObservations(..), KernelFunction, Matrix, unInputObs)
import GPLVM.Util (randomMatrixD)

type GPMean a = Matrix D a -> Matrix D a

data GaussianProcess a = GaussianProcess
    { _kernelGP :: KernelFunction a
    , _meanGP ::  GPMean a
    }

data GPTrainingData a = GPTrainingData
    { _inputTrain :: Matrix D a
    , _outputTrain :: Matrix D a
    }

makeLenses ''GaussianProcess
makeLenses ''GPTrainingData

newtype PosteriorSample a = PosteriorSample { unSample :: Matrix D a }

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
    -> Int               --- number of samples for prior
    -> Matrix D a
functionalPrior matrix@(ADelayed (Z :. rows :. cols) _) gen sampleNumber =
    delay $ matrix `mulS` randomCoeffs
    where
        randomCoeffs = randomMatrixD gen (cols, sampleNumber)

gpToPosteriorSample
    :: forall a g.
    ( GPConstraint a
    , Eq a
    , RandomGen g
    )
    => InputObservations a
    -> GaussianProcess a
    -> GPTrainingData a
    -> g
    -> Int
    -> Maybe (PosteriorSample a)
gpToPosteriorSample inputObserve gP trainingData gen sampleNumber = do
    let covarianceMatrix = testValueCovariaceMatrix gP inputObserve
        -- get covariance matrix for test inputs
    let inputTrain' = trainingData ^. inputTrain

        -- kernel applied to input training points
    let trainingKernel = (gP ^. kernelGP) inputTrain' inputTrain'
        -- Cholesky decomposition applied to kernel of training points
    let cholK = cholSH trainingKernel

        -- test points mean
    let testPointMean = (gP ^. kernelGP) (inputObserve ^. unInputObs) inputTrain'

    cholKSolve <- linearSolveS cholK testPointMean

        -- solve linear system for output training points
    cholKSolveOut <- linearSolveS cholK (trainingData ^. outputTrain)

    let cholKSolve' = delay cholKSolve
    let cholKSolveOut' = delay cholKSolveOut
    let mulD m n = delay (mulS m n)

        -- compute mean
    let mean = transpose cholKSolve' `mulD` cholKSolveOut'

        -- posterior
    let postF' = cholSH $
                     covarianceMatrix -^
                     (transpose cholKSolve' `mulD` cholKSolve')
    let prior = functionalPrior postF' gen sampleNumber
    return $ PosteriorSample $ mean +^ prior
