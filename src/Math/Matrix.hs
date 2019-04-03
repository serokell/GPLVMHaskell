{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE FlexibleInstances #-}

module Math.Matrix
       ( DimVector (..)
       , DimMatrix (..)
       , GPConstraint
       , withVec
       , withMat
       , toMatrix
       , mulM
       , transposeM
       , mapMM
       , mapDiagonalM
       , invSM
       , substractMeanM
       , (+^^)
       , (-^^)
       , (*^^)
       , cholM
       , foldAllSM
       , sumAllSM
       , trace2SM
       , toDimMatrix
       , trace2SM
       , detSM
       , toMatrixM'
       , zipWithDim
       , zipWithArray
       , identM
       , delayMatrix
       , linearSolveM
       , vectorLength
       , randomMatrixD
       , matrixRowsNum
       , matrixColsNum
       , pinvSM
       , mulPM
       , eigSHM
       )
       where

import Universum hiding (Vector, transpose, map, natVal,
                         zipWith)

import           GHC.TypeLits hiding (someNatVal)

import           Data.Array.Repa
import           Data.Array.Repa.Algorithms.Matrix hiding (trace2S)
import           Data.Random.Normal (normals)
import           Data.Vector.Unboxed.Base (Unbox)
import           Data.Vinyl.TypeLevel (AllConstrained)
import           Numeric.Dimensions
import           Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import           System.Random (Random, RandomGen)

import GPLVM.Types
import GPLVM.Util hiding (randomMatrixD, zipWithArray)

newtype DimVector r (n :: Nat) a
  = DimVector { runVector :: Vector r a }
  deriving Eq

withVec
  :: Vector D a
  -> (forall n. KnownNat n => DimVector D n a -> k)
  -> k
withVec vec f =
  case someNatVal (fromIntegral x) of
    SomeNat (Proxy :: Proxy m) -> f (DimVector @D @m vec)
  where
    (Z :. x) = extent vec

newtype DimMatrix r (y :: Nat) (x :: Nat) a
  = DimMatrix { getInternal :: Matrix r a}

instance Functor (DimVector D n) where
  fmap f (DimVector vector) = DimVector (smap f vector)

withMat
  :: Matrix D a
  -> (forall x y. (KnownNat x, KnownNat y) => DimMatrix D x y a -> k)
  -> k
withMat m f =
  let (Z :. x :. y) = extent m
  in
  case someNatVal (fromIntegral x) of
    SomeNat (Proxy :: Proxy m) -> case someNatVal (fromIntegral y) of
      SomeNat (Proxy :: Proxy n) -> f (DimMatrix @_ @m @n m)

mulM
  :: forall y1 x1 y2 x2 a.
  ( AllConstrained KnownNat [x1, x2, y1, y2]
  , Numeric a
  , x1 ~ y2
  )
  => DimMatrix D y1 x1 a
  -> DimMatrix D y2 x2 a
  -> DimMatrix D y1 x2 a
mulM (DimMatrix m1) (DimMatrix m2) = DimMatrix $ delay $ m1 `mulS` m2

transposeM
  :: (KnownNat y, KnownNat x)
  => DimMatrix D y x a
  -> DimMatrix D x y a
transposeM (DimMatrix m) = DimMatrix $ transpose m

mapMM
  ::
  ( KnownNat y
  , KnownNat x
  , Unbox a
  , Unbox b
  )
  => (a -> b)
  -> DimMatrix D y x a
  -> DimMatrix D y x b
mapMM f (DimMatrix m) =  DimMatrix $ map f m

mapDiagonalM
  ::
  ( KnownNat y
  , KnownNat x
  , Unbox a
  )
  => (a -> a)
  -> DimMatrix D y x a
  -> DimMatrix D y x a
mapDiagonalM f (DimMatrix m) = DimMatrix $ mapDiagonal f m

invSM
  ::
  ( KnownNat y
  , KnownNat x
  , Field a
  , Numeric a
  , y ~ x
  )
  => DimMatrix D y x a
  -> DimMatrix D y x a
invSM (DimMatrix m) = DimMatrix $ delay $ invS m

substractMeanM
  ::
  ( KnownNat y
  , KnownNat x
  )
  => DimMatrix D y x Double
  -> DimMatrix D y x Double
substractMeanM (DimMatrix m) = DimMatrix $ substractMean m

infixl 6 +^^, -^^
infixl 7 *^^

(+^^)
  :: forall y1 x1 y2 x2 a.
  ( AllConstrained KnownNat [x1, x2, y1, y2]
  , x1 ~ x2
  , y1 ~ y2
  , Num a
  )
  => DimMatrix D y1 x1 a
  -> DimMatrix D y2 x2 a
  -> DimMatrix D y2 x2 a
(+^^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 +^ m2

(-^^)
  :: forall y1 x1 y2 x2 a.
  ( AllConstrained KnownNat [x1, x2, y1, y2]
  , x1 ~ x2
  , y1 ~ y2
  , Num a
  )
  => DimMatrix D y1 x1 a
  -> DimMatrix D y2 x2 a
  -> DimMatrix D y2 x2 a
(-^^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 -^ m2

(*^^)
  :: forall y1 x1 y2 x2 a.
  ( AllConstrained KnownNat [x1, x2, y1, y2]
  , x1 ~ x2
  , y1 ~ y2
  , Num a
  )
  => DimMatrix D y1 x1 a
  -> DimMatrix D y2 x2 a
  -> DimMatrix D y2 x2 a
(*^^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 *^ m2

cholM
  ::
  ( KnownNat y
  , KnownNat x
  , Field a
  , y ~ x
  )
  => DimMatrix D y x a
  -> DimMatrix D y x a
cholM (DimMatrix m) = DimMatrix $ delay $ chol $ trustSym $ computeS m

sumAllSM
  ::
  ( KnownNat y
  , KnownNat x
  , Num a)
  => DimMatrix D y x a
  -> a
sumAllSM (DimMatrix m) = sumAllS m

foldAllSM
  :: (KnownNat y, KnownNat x)
  => (Double -> Double -> Double)
  -> Double
  -> DimMatrix D y x Double
  -> Double
foldAllSM f initValue (DimMatrix m) = foldAllS f initValue m

detSM
  :: (KnownNat y, KnownNat x)
  => DimMatrix D y x Double
  -> Double
detSM (DimMatrix m) = detS m

trace2SM
  :: (KnownNat x, KnownNat y)
  => DimMatrix D x y Double
  -> Double
trace2SM (DimMatrix m) = trace2S $ computeS m

toDimMatrix
  ::
  ( Source r a
  , KnownNat m
  , KnownNat n
  )
  => DimVector r m a
  -> Int
  -> DimMatrix D m n a
toDimMatrix (DimVector arr) desiredSize =
  DimMatrix (toMatrix arr desiredSize)

toMatrixM'
  :: (Source r a, KnownNat n, KnownNat m)
  => DimVector r n a
  -> DimMatrix D n m a
toMatrixM' (DimVector arr) =
  DimMatrix $ fromFunction (Z :. dimension :. 1) generator
  where
    dimension = size . extent $ arr
    generator (Z :. rows :. cols) = linearIndex arr rows

zipWithDim
  ::
  ( Source r1 a
  , Source r2 b
  , KnownNat m
  , KnownNat n
  )
  => (a -> b -> c)
  -> DimMatrix r1 m n a
  -> DimMatrix r2 m n b
  -> DimMatrix D m n c
zipWithDim f (DimMatrix mat1) (DimMatrix mat2) =
  DimMatrix $ zipWith f mat1 mat2

zipWithArray
  ::
  ( KnownNat n
  , KnownNat m
  )
  => (a -> b -> c)
  -> DimVector D n a
  -> DimMatrix D n m b
  -> DimMatrix D n m c
zipWithArray f (DimVector array1) (DimMatrix array2) =
  DimMatrix $ zipWith f (toMatrix array1 n) array2
  where
    (Z :. n :. _) = extent array2

identM
  :: forall m n a.
  ( KnownNat m
  , KnownNat n
  , GPConstraint a
  , m ~ n
  )
  => DimMatrix D m n a
identM =
  let dim = fromEnum $
            natVal @m @Proxy Proxy in DimMatrix $
            identD dim

delayMatrix
  ::
  ( KnownNat m
  , KnownNat n
  , Source r a
  )
  => DimMatrix r m n a
  -> DimMatrix D m n a
delayMatrix (DimMatrix matrix) = DimMatrix (delay matrix)

linearSolveM
  ::
  ( Field a
  , AllConstrained KnownNat '[m, n]
  )
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

matrixRowsNum
  :: forall r m n a. KnownNat m
  => DimMatrix r m n a
  -> Int
matrixRowsNum _ = fromEnum $ natVal (Proxy @m)

matrixColsNum
  :: forall r m n a. KnownNat n
  => DimMatrix r m n a
  -> Int
matrixColsNum _ = fromEnum $ natVal (Proxy @n)

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
    DimMatrix . delay $ fromListUnboxed (Z :. rows :. cols) randomList

type GPConstraint a =
  ( Field a
  , Random a
  , Unbox a
  , Floating a
  , Eq a
  )

hermToMatrixM
  :: Herm a
  -> DimMatrix D m n a
hermToMatrixM = undefined

pinvSM
  ::
  ( KnownNat y
  , KnownNat x
  , Field a
  , Numeric a
  , y ~ x
  )
  => DimMatrix D y x a
  -> DimMatrix D y x a
pinvSM (DimMatrix m) = DimMatrix . delay $ pinvS m

mulPM
  :: forall y1 x1 y2 x2 a.
  ( AllConstrained KnownNat [x1, x2, y1, y2]
  , Numeric a
  , x1 ~ y2
  )
  => DimMatrix D y1 x1 a
  -> DimMatrix D y2 x2 a
  -> DimMatrix D y1 x2 a
mulPM (DimMatrix m) (DimMatrix n) =
  DimMatrix . delay . runIdentity $ m `mulP` n

eigSHM
  :: (Field a, Numeric a)
  => Herm a
  -> (DimVector D m Double, DimMatrix D m m a)
eigSHM hermM =
  (DimVector $ delay $ fst out, DimMatrix $ delay $ snd out)
  where
    out = eigSH hermM
