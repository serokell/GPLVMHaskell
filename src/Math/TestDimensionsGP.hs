{-# LANGUAGE ConstraintKinds #-}

module Math.TestDimensionsGP where

import Universum hiding (Vector)

import Control.Lens (makeLenses)

import Data.Array.Repa (Array (..), D, DIM1, DIM2, Source, Z (..), (:.) (..))
import qualified Data.Array.Repa as R
import GHC.TypeLits hiding (someNatVal)
import GPLVM.Types
import Numeric.Dimensions

import Data.Array.Repa.Repr.Unboxed (Unbox)
import Data.Array.Repa.Algorithms.Matrix
import Data.Vinyl.TypeLevel (AllConstrained)

import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

import Math.TestDimensions

import System.Random


newtype DimVector r (x :: Nat) a = DimVector { runVector :: Vector r a }

withVec
  :: Vector D a
  -> (forall n. KnownNat n => DimVector D n a -> k)
  -> k
withVec vec f =
  case someNatVal (fromIntegral x) of
    SomeNat (Proxy :: Proxy m) -> f (DimVector @D @m vec)
  where
    (Z :. x) = R.extent vec
{-
toMatrix
  :: Source r a
  => Array r DIM1 a
  -> Int
  -> Array D DIM2 a
toMatrix arr desiredSize =
  R.fromFunction (Z :. desiredSize :. newDimension) generator
  where
    generator (Z :. rows :. cols) = R.linearIndex arr cols
    newDimension = R.size $ R.extent $ arr

toDimMatrix
  :: Source r a
  => DimVector r m a
  -> Int
  -> DimMatrix D m n a
toDimMatrix (DimVector arr) desiredSize =
  DimMatrix (toMatrix arr desiredSize)

toMatrix'
  :: (Source r a, KnownNat n, KnownNat m)
  => DimVector r n a
  -> DimMatrix D n m a
toMatrix' (DimVector arr) =
  DimMatrix $ R.fromFunction (Z :. dimension :. 1) generator
  where
    dimension = R.size . R.extent $ arr
    generator (Z :. rows :. cols) = R.linearIndex arr rows

zipWithDim
  :: (Source r1 a, Source r2 b, KnownNat m, KnownNat n)
  => (a -> b -> c)
  -> DimMatrix r1 m n a
  -> DimMatrix r2 m n b
  -> DimMatrix D m n c
zipWithDim f (DimMatrix mat1) (DimMatrix mat2) =
  DimMatrix $ R.zipWith f mat1 mat2

zipWithArray
  :: (KnownNat n, KnownNat m)
  => (a -> b -> c)
  -> DimVector D n a
  -> DimMatrix D n m b
  -> DimMatrix D n m c
zipWithArray f (DimVector array1) (DimMatrix array2) =
  DimMatrix $ R.zipWith f (toMatrix array1 n) array2
  where
     (Z :. n :. _) = R.extent array2

zipWithArray'
  :: (a -> b -> c)
  -> Array D DIM2 a
  -> Array D DIM1 b
  -> Array D DIM2 c
zipWithArray' f array1@(ADelayed (Z :. (n :: Int) :. _) _) array2 =
  R.zipWith f array1 (toMatrix array2 n)

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


newtype GaussianProcess a = KernelFunction
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

-- | Main GP function: get a posterior sample by some kernel function, input observations and training data

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
-}
