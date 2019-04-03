{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleInstances #-}

module Math.TestDimensionsGP where

import Universum hiding (Vector)

import Control.Lens (makeLenses)

import Data.Array.Repa ((:.)(..), Array(..), D, DIM1, DIM2, Source, Z(..))
import qualified Data.Array.Repa as R
import GHC.TypeLits hiding (someNatVal)
import GPLVM.Types hiding (InputObservations)
import Numeric.Dimensions

import Data.Array.Repa.Algorithms.Matrix
import Data.Array.Repa.Repr.Unboxed (Unbox)
import Data.Vinyl.TypeLevel (AllConstrained)

import Data.Random.Normal (normals)
import Data.Vinyl.TypeLevel (AllConstrained)

--import Math.TestDimensions

import System.Random

import Math.Matrix

newtype InputObservations (n :: Nat) a
  = InputObservations { _unInputObs :: DimVector D n a }
  deriving Eq

makeLenses ''InputObservations

{-
infixl 6 ++^
infixl 6 --^
infixl 7 **^
infixl 7 //^

(++^), (**^), (--^), (//^)
  :: (Num a, Fractional a)
  => Array D DIM1 a
  -> Array D DIM2 a
  -> Array D DIM2 a
(++^) = zipWithArray (+)
(**^) = zipWithArray (*)
(--^) = zipWithArray (-)
(//^) = zipWithArray (/)

infixl 6 ^++
infixl 6 ^--
infixl 7 ^**
infixl 7 ^//


(^++), (^**), (^--), (^//)
   :: (Num a, Fractional a)
   => Array D DIM2 a
   -> Array D DIM1 a
   -> Array D DIM2 a
(^++) = zipWithArray' (+)
(^**) = zipWithArray' (*)
(^--) = zipWithArray' (-)
(^//) = zipWithArray' (/)
-}

newtype GaussianProcess (n :: Nat) (m :: Nat) a = GaussianProcess
    { -- | we define gaussian process by some kernel function
      _kernelGP :: DimVector D n a -> DimVector D m a -> DimMatrix D n m a
    }

data GPTrainingData (n :: Nat) a = GPTrainingData
    { _inputTrain  :: DimVector D n a  -- ^ input training data
    , _outputTrain :: DimVector D n a  -- ^ output training data
    }

makeLenses ''GaussianProcess
makeLenses ''GPTrainingData

-- | Main GP function: get a posterior sample by some kernel function, input observations and training data

gpToPosteriorSample
  :: forall a m p.
  ( GPConstraint a
  , AllConstrained KnownNat [m, p]
  )
  => InputObservations m a        -- ^ input observations
  -> (forall c d. DimVector D c a -> DimVector D d a -> DimMatrix D c d a)  -- ^ kernel function
  -> GPTrainingData p a           -- ^ training data
  -> Int                          -- ^ number of samples
  -> Maybe (Matrix D a)           -- ^ posterior functional prior
gpToPosteriorSample (InputObservations observe) kernel trainingData sampleNumber = do
  -- | kernel applied to input test points (so-called K_ss)
  let covarianceMatrix = kernel @m @m observe observe

  -- | kernel applied to input training points (so-called K)
  let trainingKernel = kernel @p @p inputTrain' inputTrain'

  -- | Cholesky decomposition applied to kernel of training points (:), so-called L
  let cholK = cholM $ trainingKernel +^^ (mapMM (* 0.00005) (identM @p))

  -- | covariance between test points and input training points (so-called K_s)
  let testPointMean = kernel @p @m inputTrain' observe

  -- | (roots of L * x = K_s)
  cholKSolve <- linearSolveM cholK testPointMean

  -- | solve linear system for output training points
  cholKSolveOut <- linearSolveM cholK (transposeM $ toDimMatrix outputTrain' len)

  -- | compute mean
  let mean = (transposeM cholKSolve) `mulM` cholKSolveOut

  -- | posterior
  let postF = cholM $
                     covarianceMatrix +^^
                     mapMM (* 1.0e-6) (identM @m) -^^
                     (transposeM cholKSolve) `mulM` cholKSolve
  -- | posterior sample
  return $ (getInternal $ mean +^^ functionalPrior postF)
    where
      inputTrain' = trainingData ^. inputTrain
      outputTrain' = trainingData ^. outputTrain
      len = vectorLength observe
      functionalPrior matrix =
        matrix `mulM` randomCoeffs
        where
          randomCoeffs = randomMatrixD (mkStdGen 0) (rows, sampleNumber)
          rows = matrixRowsNum matrix
