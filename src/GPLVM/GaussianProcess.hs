{-# LANGUAGE AllowAmbiguousTypes, ConstraintKinds #-}

module GPLVM.GaussianProcess
       ( GaussianProcess (..)
       , GPTrainingData (..)
       , PosteriorSample (..)
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
    { -- | we define gaussian process by some kernel function
      _kernelGP :: Vector D a -> Vector D a -> Matrix D a
    }

data GPTrainingData a = GPTrainingData
    { _inputTrain  :: Vector D a  -- ^ input training data
    , _outputTrain :: Vector D a  -- ^ output training data
    }

makeLenses ''GaussianProcess
makeLenses ''GPTrainingData

newtype PosteriorSample a = PosteriorSample
    { unSample :: Matrix D a  -- ^ posterior sample
    }

-- | Constraint kind required to getting a posterior sample by a given GP
type GPConstraint a =
  ( Field a
  , Random a
  , Unbox a
  , Floating a
  , Eq a
  )

gpToPosteriorSample
  :: forall a g.
  GPConstraint a
  => InputObservations a        -- ^ input observations
  -> GaussianProcess a          -- ^ kernel function
  -> GPTrainingData a           -- ^ training data
  -> Int                        -- ^ number of samples
  -> Maybe (PosteriorSample a)  -- ^ posterior functional prior
gpToPosteriorSample (InputObservations observe@(ADelayed (Z :. len) _)) gP trainingData sampleNumber = do
  -- | kernel applied to input test points (so-called K_ss)
  let covarianceMatrix = kernel observe observe

  -- | kernel applied to input training points (so-called K)
  let trainingKernel = kernel inputTrain' inputTrain'

  -- | Cholesky decomposition applied to kernel of training points (:), so-called L
  let cholK = cholSH $ trainingKernel +^
                (smap (* 0.00005) . identD . size . extent $ inputTrain')

  -- | covariance between test points and input training points (so-called K_s)
  let testPointMean = kernel observe inputTrain'

  -- | (roots of L * x = K_s)
  cholKSolve <- delay <$> linearSolveS cholK testPointMean

  -- | solve linear system for output training points
  cholKSolveOut <- delay <$> linearSolveS cholK (transposeMatrix $ toMatrix outputTrain' len)

  -- | compute mean
  let mean = (transposeMatrix cholKSolveOut) `mulD` cholKSolve

  -- | posterior
  let postF' = cholSH $
                     covarianceMatrix +^
                     ((smap (* 1.0e-6) (identD len)) -^
                     (transposeMatrix cholKSolve) `mulD` cholKSolve)

  -- | posterior sample
  return $ (PosteriorSample $ mean +^ (functionalPrior postF' (mkStdGen 1) sampleNumber))
    where
       kernel = gP ^. kernelGP
       inputTrain' = trainingData ^. inputTrain
       outputTrain' = trainingData ^. outputTrain
       mulD m n = delay $ m `mulS` n
       functionalPrior matrix@(ADelayed (Z :. rows :. cols) _) gen sampleNumber =
         delay $ matrix `mulS` randomCoeffs
         where
           randomCoeffs = randomMatrixD (mkStdGen (-4)) (rows, sampleNumber)
