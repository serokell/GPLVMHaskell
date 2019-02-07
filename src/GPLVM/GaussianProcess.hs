{-# LANGUAGE AllowAmbiguousTypes, ConstraintKinds #-}

module GPLVM.GaussianProcess
       ( GaussianProcess (..)
       , GPTrainingData (..)
       , PosteriorSample (..)
       , functionalPrior
       , gpToPosteriorSample
       , kernelGP
       ) where

import Universum hiding (transpose, Vector)

import Control.Lens (makeLenses)

import Data.Array.Repa
import Data.Array.Repa.Repr.Unboxed (Unbox)
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import System.Random (Random, RandomGen, mkStdGen)

import GPLVM.Types (InputObservations(..), KernelFunction, Matrix,
                    Vector, unInputObs)
import GPLVM.Util

data GaussianProcess a = GaussianProcess
    { _kernelGP :: Vector D a -> Vector D a -> Matrix D a
    }

data GPTrainingData a = GPTrainingData
    { _inputTrain :: Vector D a
    , _outputTrain :: Vector D a
    }

makeLenses ''GaussianProcess
makeLenses ''GPTrainingData

newtype PosteriorSample a = PosteriorSample { unSample :: Matrix D a }

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
        randomCoeffs = randomMatrixD gen (rows, sampleNumber)

gpToPosteriorSample
    :: forall a g.
    ( GPConstraint a
    , Eq a
    )
    => InputObservations a
    -> GaussianProcess a
    -> GPTrainingData a
    -> Int
    -> Maybe (PosteriorSample a, Vector D a)  -- returns posterior functional prior and standard deviation
gpToPosteriorSample (InputObservations observe@(ADelayed (Z :. len) _)) gP trainingData sampleNumber = do
        -- kernel applied to input test points (so-called K_ss)
    let covarianceMatrix = kernel observe observe

        -- kernel applied to input training points (so-called K)
    let trainingKernel = kernel inputTrain' inputTrain'

        -- Cholesky decomposition applied to kernel of training points (:), so-called L
    let cholK = cholSH (trainingKernel +^ ((delay . smap (* 0.00005)) $
                        ident . size . extent $ inputTrain'))

        -- covariance between test points and input training points (so-called K_s)
    let testPointMean = kernel observe inputTrain'

        -- (roots of L * x = K_s)
    cholKSolve <- delay <$> linearSolveS cholK testPointMean

        -- solve linear system for output training points
    cholKSolveOut <- transposeMatrix . delay <$> linearSolveS cholK (toMatrix' $ outputTrain')
        -- compute mean
    let mean = cholKSolveOut `mulD` cholKSolve

        -- compute standard deviation

    let standardDeviation = ((delay . smap sqrt) $ (takeDiagD covarianceMatrix)) -^ ((delay . sumS) $ smap flipExp (transposeMatrix cholKSolve))

        -- posterior
    let postF' = cholSH $
                     covarianceMatrix +^ ((smap (* 1.0e-6) (identD len)) -^
                                          (transposeMatrix cholKSolve) `mulS` cholKSolve)
    let prior = functionalPrior postF' (mkStdGen 0) sampleNumber
    return $ (PosteriorSample $ mean +^ prior, standardDeviation)
    where
       kernel = gP ^. kernelGP
       inputTrain' = trainingData ^. inputTrain
       outputTrain' = trainingData ^. outputTrain
       mulD m n = delay $ m `mulS` n
       flipExp = (flip (^)) 2
