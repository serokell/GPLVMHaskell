{-# LANGUAGE AllowAmbiguousTypes #-}

module GPLVM.TypeSafe.GPLVM
  ( GPLVMTypeSafe (..)
  , GPLVM (..)
  , WithIVM (..)
  , makeGPLVM
  , makeGPLVMTypeSafe
  ) where

import Universum hiding (Nat, (%~))

import           Data.Array.Repa hiding (Z)
import           Data.Default (Default (def))
import           Data.Singletons.Decide (Decision (..), (:~:) (..), (%~))
import           Data.Type.Natural
import           Data.Vinyl.TypeLevel (AllConstrained)

import           GPLVM.TypeSafe.Util (withMat, Decisions)
import           GPLVM.TypeSafe.Types (DimMatrix (..), (:<:) (..), (%<))
import           GPLVM.Types (Matrix)

import           Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

data GPLVMTypeSafe =
  forall (d :: Nat)
  (m :: Nat)
  (n :: Nat).
  (d <= n ~ 'True) => GPLVMTypeSafe
  { iterations  :: Int
  , covGPLVM    :: DimMatrix D n n Double
  , inputGPLVM  :: DimMatrix D m n Double
  , resultGPLVM :: DimMatrix D m d Double
  }

data GPLVM = GPLVM
  { iter         :: Int
  , covarUnsafe  :: Matrix D Double
  , inputUnsafe  :: Matrix D Double
  , resultUnsafe :: Matrix D Double
  }

data WithIVM =
    Yes
  | No
  deriving (Show, Read)

instance Default WithIVM where
  def = No

makeGPLVMTypeSafe
  :: forall (d :: Nat) (m1 :: Nat) (m2 :: Nat) (n1 :: Nat) (n2 :: Nat).
  ( AllConstrained SingI [m1, m2, n1, n2, d]
  , m1 ~ n1
  , n1 ~ n2
  )
  => Int                       -- ^ the number of iterations
  -> DimMatrix D m1 n1 Double  -- ^ covariance
  -> DimMatrix D m2 n2 Double  -- ^ input
  -> WithIVM
  -> GPLVMTypeSafe
makeGPLVMTypeSafe = undefined

makeGPLVM
  :: Int              -- ^ a number of active points
  -> Int              -- ^ a number of iterations
  -> Matrix D Double  -- ^ covariance matrix
  -> Matrix D Double  -- ^ input matrix
  -> WithIVM
  -> GPLVM
makeGPLVM active iterations cov input ivmFlag =
  case toSing (intToNat active) of
    SomeSing (active :: Sing active) -> withSingI active $ withMat cov $
      \(cov' :: DimMatrix D y x Double) -> withMat input $
      \(input' :: DimMatrix D y1 x1 Double) ->
      let (sol1, sol2, sol3) = checkConditions @y @x @y1 @x1 @active in
      case (sol1, sol2) of
        _ -> error "equalities are false"
        (Proved Refl, Proved Refl) -> case sol3 of
          Proved LS -> fromTypeSafeToGPLVM $ makeGPLVMTypeSafe @active iterations cov' input' ivmFlag
          _ -> error "unequality is false"
  where
    checkConditions
      :: forall (y :: Nat) (x :: Nat) (y2 :: Nat) (x2 :: Nat) (d :: Nat).
      ( AllConstrained SingI '[x, y, x2, d]
      )
      => Decisions (x :~: y, x :~: x2, d :<: x)
    checkConditions =
      let des1 = (sing :: Sing x) %~ (sing :: Sing y)
          des2 = (sing :: Sing x) %~ (sing :: Sing x2)
          des3 = (sing :: Sing d) %< (sing :: Sing x)
      in (des1, des2, des3)


    fromTypeSafeToGPLVM
      :: GPLVMTypeSafe
      -> GPLVM
    fromTypeSafeToGPLVM (GPLVMTypeSafe iter cov input result) =
      GPLVM iter cov' input' result'
      where
        cov' = getInternal cov
        input' = getInternal input
        result' = getInternal result
