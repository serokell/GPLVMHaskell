{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleInstances #-}

module Math.TestDimensionsGP where

import Universum hiding (Vector)

import Control.Lens (makeLenses)

import Data.Array.Repa (Array (..), D, DIM1, DIM2, Source, Z (..), (:.) (..))
import qualified Data.Array.Repa as R
import GHC.TypeLits hiding (someNatVal, natVal)
import GPLVM.Types hiding (InputObservations)
import Numeric.Dimensions

import Data.Array.Repa.Repr.Unboxed (Unbox)
import Data.Array.Repa.Algorithms.Matrix

import Data.Random.Normal (normals)
import Data.Vinyl.TypeLevel (AllConstrained)

import qualified Data.Vector.Storable as V

import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

import Math.TestDimensions

import System.Random


newtype DimVector r (n :: Nat) a
  = DimVector { runVector :: Vector r a }
  deriving Eq

newtype InputObservations (n :: Nat) a
  = InputObservations { _unInputObs :: DimVector D n a }
  deriving Eq

makeLenses ''InputObservations

withVec
  :: Vector D a
  -> (forall n. KnownNat n => DimVector D n a -> k)
  -> k
withVec vec f =
  case someNatVal (fromIntegral x) of
    SomeNat (Proxy :: Proxy m) -> f (DimVector @D @m vec)
  where
    (Z :. x) = R.extent vec

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
  :: (Source r a, KnownNat m, KnownNat n)
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

instance Functor (DimVector D n) where
    fmap f (DimVector vector) = DimVector (R.smap f vector)

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

newtype PosteriorSample a = PosteriorSample
    { unSample :: forall m n . DimMatrix D m n a  -- ^ posterior sample
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

identM
  :: (GPConstraint a, KnownNat m, KnownNat n, m ~ n)
  => Int
  -> DimMatrix D m n a
identM dim = DimMatrix $ identD dim

delayMatrix
  :: (KnownNat m, KnownNat n, Source r a)
  => DimMatrix r m n a
  -> DimMatrix D m n a
delayMatrix (DimMatrix matrix) = DimMatrix (R.delay matrix)

linearSolveM
  :: (Field a, KnownNat m, KnownNat n)
  => DimMatrix D m m a
  -> DimMatrix D m n a
  -> Maybe (DimMatrix D m n a)
linearSolveM (DimMatrix mat1) (DimMatrix mat2) =
  case linearSolveS mat1 mat2 of
    Nothing  -> Nothing
    Just sol -> Just . delayMatrix . DimMatrix $ sol

vectorLength
  :: forall r m a. KnownNat m
  => DimVector r m a
  -> Int
vectorLength _ = fromEnum $ natVal (Proxy @m)

matrixDims
  :: forall r m n a.
  ( KnownNat m
  , KnownNat n
  )
  => DimMatrix r m n a
  -> (Int, Int)
matrixDims _ = fromEnumPair (Proxy @m, Proxy @n)
  where
    fromEnumPair (x, y) = (fromEnum x, fromEnum y)


randomMatrixD
    :: forall a g m n.
    ( RandomGen g
    , Random a
    , Unbox a
    , Floating a
    , KnownNat m
    , KnownNat n
    )
    => g
    -> (Int, Int)
    -> DimMatrix D m n a
randomMatrixD gen (rows, cols) =
    let randomList = take (rows * cols) (normals gen) in
    DimMatrix . R.delay $ R.fromListUnboxed (Z :. rows :. cols) randomList

gpToPosteriorSample
  :: forall a m n o p.
  ( GPConstraint a
  , AllConstrained KnownNat [m, n, o, p]
  , o ~ n
  , p ~ n
  , m ~ n
  )
  => InputObservations m a        -- ^ input observations
  -> GaussianProcess n o a        -- ^ kernel function
  -> GPTrainingData p a           -- ^ training data
  -> Int                          -- ^ number of samples
  -> Maybe (PosteriorSample a)    -- ^ posterior functional prior
gpToPosteriorSample (InputObservations observe) gP trainingData sampleNumber = do
  -- | kernel applied to input test points (so-called K_ss)
  let covarianceMatrix = kernel observe observe

  -- | kernel applied to input training points (so-called K)
  let trainingKernel = kernel inputTrain' inputTrain'

  -- | Cholesky decomposition applied to kernel of training points (:), so-called L
  let cholK = cholM $ trainingKernel +^^
                (mapMM (* 0.00005) . identM $ vectorLength inputTrain')

  -- | covariance between test points and input training points (so-called K_s)
  let testPointMean = kernel observe inputTrain'

  -- | (roots of L * x = K_s)
  cholKSolve <- linearSolveM cholK testPointMean

  -- | solve linear system for output training points
  cholKSolveOut <- linearSolveM cholK (transposeM $ toDimMatrix outputTrain' len)

  -- | compute mean
  let mean = (transposeM cholKSolveOut) *^^ cholKSolve

  -- | posterior
  let postF' = cholM $
                     covarianceMatrix +^^
                     ((mapMM (* 1.0e-6) (identM len)) -^^
                     (transposeM cholKSolve) *^^ cholKSolve)

  -- | posterior sample
  return $ (PosteriorSample $ mean +^^ (functionalPrior postF' sampleNumber))
    where
       kernel = gP ^. kernelGP
       inputTrain' = trainingData ^. inputTrain
       outputTrain' = trainingData ^. outputTrain
       len = vectorLength observe
       functionalPrior matrix sampleNumber =
         matrix *^^ randomCoeffs
         where
           randomCoeffs = randomMatrixD (mkStdGen (-4)) (rows, sampleNumber)
           rows = fst $ matrixDims matrix
