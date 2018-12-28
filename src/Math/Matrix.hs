{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE PolyKinds #-}

module Math.Matrix where

import Universum hiding (someNatVal, natVal)

import Data.Array.Repa hiding (rank)
import Data.Array.Repa.Repr.ForeignPtr
import Data.Bifunctor (first)
import Data.Vinyl.TypeLevel (AllConstrained, Fst, Snd)
import Data.Singletons
import Data.Singletons.Sigma
import Data.Singletons.TH
import Data.Proxy
import Numeric.LinearAlgebra.HMatrix (Complex, Field, LSDiv,
                                      Normed, Numeric, Product,
                                      RandDist(..), RealElement, Seed,
                                      Vector)
import qualified Numeric.LinearAlgebra.HMatrix as H
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import Unsafe.Coerce

import GHC.TypeLits

data Matrix r (m :: Nat) (n :: Nat) :: Type -> Type where
    Matrix :: Array r DIM2 a -> Matrix r m n a

unMatrix :: Matrix r m n a -> Array r DIM2 a
unMatrix (Matrix matrix) = matrix

type KnownNats xs = AllConstrained KnownNat xs

data SomeMatrix r a where
   SomeMatrix :: KnownNats [m, n] => Matrix r m n a -> SomeMatrix r a

type MatrixD = Matrix D

--type family (a,b) :: Type

type ExistsMatrixD a = Sigma (Nat,Nat)

{-
makeMatrix ::
    forall a.
    Array D DIM2 a ->
    exists n m.
    (KnownNat n, KnownNat m, Matrix n m a)
makeMatrix = -}

makeMatrix
    :: forall r a. Source r a
    => Array r DIM2 a
    -> SomeMatrix r a
makeMatrix array =
  let (Z :. y :. x) = extent array
      (Just k) = someNatVal $ fromIntegral y in
      case k of
      	SomeNat (Proxy :: Proxy n) -> case (someNatVal $ fromIntegral y) of
      	  Just (SomeNat (Proxy :: Proxy m)) -> SomeMatrix $ (Matrix array :: Matrix r m n a)
          Nothing -> error "izviniti: pepe_sad"


mulM
    :: forall x1 x2 y1 y2 .
    ( KnownNats [x1, x2, y1, y2]
    , y1 ~ x2
    )
    => Matrix D x1 y1 Double
    -> Matrix D x2 y2 Double
    -> Matrix F x1 y2 Double
mulM mk1 mk2 =
  let res@(AForeignPtr (Z :. y :. x) _ _) = (unMatrix mk1) `mulS` (unMatrix mk2)
  in if (natVal (Proxy @x1) == fromIntegral x) && (natVal (Proxy @y1) == fromIntegral y)
  	 then Matrix $ res
  	 else (error "test")

{-
sMult
    :: forall m2 n2 m3 n3 .
    KnownNats [m2, n2, m3, n3]
    => SomeMatrix D Double
    -> SomeMatrix D Double
    -> SomeMatrix F Double
sMult matr1 matr2 =
	case matr1 of
		(SomeMatrix (mat1 :: forall m1 n1 . (KnownNat n1, KnownNat m1) => Matrix D m1 n1 Double)) -> case matr2 of
		  (SomeMatrix (mat2 :: Matrix D m2 n2 Double)) -> SomeMatrix $ ( mat1 `mulM` mat2 ::  Matrix F m3 n3 Double)
-}
