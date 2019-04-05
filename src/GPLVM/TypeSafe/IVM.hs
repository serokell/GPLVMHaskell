-- | based on Informative vector machine implemented as follows:
-- https://github.com/lawrennd/ivm/blob/master/matlab/ivmCreate.m
-- See also
-- https://papers.nips.cc/paper/2240-fast-sparse-gaussian-process-methods-the-informative-vector-machine.pdf

{-# LANGUAGE AllowAmbiguousTypes #-}

module GPLVM.TypeSafe.IVM
  ( IVM (..)
  , makeIVM
  , makeIVMTypeSafe
  ) where

import Universum hiding (All, Nat, One, Vector, transpose, (%~))

import Prelude (log)

import Control.Lens (makeLenses)

import Data.Array.Repa hiding (Z)
import Data.Singletons
import Data.Singletons.Decide
import Data.Type.Natural
import Data.Vinyl.TypeLevel (AllConstrained)

import GPLVM.Types
import GPLVM.TypeSafe.Types
import GPLVM.TypeSafe.Util
import GPLVM.Util (normalDistributionProbability)


data IVMTypeSafe =
  forall (d :: Nat)                            -- d a number of active points
  (m :: Nat)                                   -- input columns
  (n :: Nat).                                  -- input rows
  (d <= n ~ 'True) => IVMInput
  { covariance    :: DimMatrix D n n Double    -- input covariance
  , inputMatrix   :: DimMatrix D m n Double    -- input matrix
  , sparsedMatrix :: DimMatrix D m d Double    -- output sparsed matrix
  , resultIndeces :: [Int]                     -- the list of row indeces in original
                                               -- | that were placed to the output sparsed matrix
  }

-- dimesionless version of IVM
data IVM = IVM
  { covar        :: Matrix D Double
  , input        :: Matrix D Double
  , sparsed      :: Matrix D Double
  , selectedRows :: [Int]
  }

type family DecisionThree k :: Type

type instance DecisionThree (a, b, c) = (Decision a, Decision b, Decision c)

data Purpose = Regression | Classification


makeIVMTypeSafe
  :: forall (d :: Nat) (m1 :: Nat) (m2 :: Nat) (n1 :: Nat) (n2 :: Nat).
  ( AllConstrained SingI [m1, m2, n1, n2, d]
  , d <= n1 ~ 'True
  , m1 ~ n1
  , n1 ~ n2
  )
  => DimMatrix D m1 n1 Double  -- covariance
  -> DimMatrix D m2 n2 Double  -- input
  -> Purpose
  -> IVMTypeSafe
makeIVMTypeSafe =
	let gIN i n = (cIN (i - 1) n) +
	      (normalDistributionProbability 0 1 (uIN (i - 1) n)) +
	      (fi (uIN (i - 1) n))
	    uIN i n = 1.0
	    cIN i n = (yN n)*(sqrt (ksiIN i n))
	    ksiIN i n = undefined
	    uIN i n = (cIN i n)*((muIN i n) + b)
	    fi = undefined
	--    normal
	in undefined

makeIVM
  :: Int              -- a number of active point
  -> Matrix D Double  -- covariance matrix
  -> Matrix D Double  -- input matrix
  -> Purpose
  -> IVM
makeIVM actPoints cov input purpose =
  case toSing (intToNat actPoints) of
    SomeSing (active :: Sing active) -> withSingI active $ withMat cov $
      \(mat1 :: DimMatrix D y x Double) -> withMat input $
      \(mat2 :: DimMatrix D y2 x2 Double) ->
      let (sol1, sol2, sol3) = checkInput @y @x @y2 @x2 @active in
      case (sol1, sol2) of
        _ -> error "equalities are false"
        (Proved Refl, Proved Refl) -> case sol3 of
          Disproved _ -> error "desired dimension is greater than required"
          Proved LS    -> convertTypeSafeIVM $ makeIVMTypeSafe @active mat1 mat2 purpose
  where
    checkInput
      :: forall (y :: Nat) (x :: Nat) (y2 :: Nat) (x2 :: Nat) (active :: Nat).
      (AllConstrained SingI '[x, y, x2, active])
      => DecisionThree (x :~: y, x :~: x2, active :<: x)
    checkInput =
      let des1 = (sing :: Sing x) %~ (sing :: Sing y)
          des2 = (sing :: Sing x) %~ (sing :: Sing x2)
          des3 = (sing :: Sing active) %< (sing :: Sing x)
      in (des1, des2, des3)

    convertTypeSafeIVM
      :: IVMTypeSafe
      -> IVM
    convertTypeSafeIVM (IVMInput (DimMatrix cov) (DimMatrix input) (DimMatrix sparsed) ind) =
      IVM cov input sparsed ind


