{-# LANGUAGE MagicHash, PolyKinds, UndecidableInstances #-}

-- | Module with main types and type families
module GPLVM.TypeSafe.Types where

import Universum hiding (Nat, One, (%~), (++))

import GHC.Exts (unsafeCoerce#)

import GPLVM.Types

import Data.Array.Repa hiding (Z)

import Data.Singletons.Decide (Decision(..))
import Data.Type.Natural

-- | Matrix data type with dimensions
newtype DimMatrix r (y :: Nat) (x :: Nat) a
  = DimMatrix { getInternal :: Matrix r a}

data A = A

-- | input and output data structure for PPCA
-- The last element in case of data without missed values should be Nothing
data PPCA = PPCA
  {  _noMissedData        :: Bool
   , _learningData        :: Matrix D Double
   , desiredDimensions    :: Int
   , stopParameter        :: Either Int Double
   , _variance            :: Double
   , _W                   :: Matrix D Double
   , _finalExpLikelihood  :: Double
   , _restoredMatrix      :: Maybe (Matrix D Double)
   }

-- | Type with unknown natural number which is less then max
data LessThen (max :: Nat) = forall (i :: Nat). (((i <= max) ~ 'True), SingI i) => Less (Proxy i)

-- | Inductively defined `+` type function for natural numbers
-- We use it because with generated by TH function from "type-naturals" library
-- GHC can't deduce some properties like (n + ('S 'Z) ~ 'S n)
type family (n :: Nat) :+: (m :: Nat) where
  (:+:) n ('Z) = n
  (:+:) n ('S l) = ('S n) :+: l


-- | Data type with evidence that a <= b
data (a :: k) :<: (b :: k) where
  LS :: forall k (a :: k) (b :: k). ((a <= b) ~ 'True) => a :<: b

-- | Class for kinds for which you can proof that one less or equal then another
class LSDecide k where
  (%<) :: forall (a :: k) (b :: k). Sing a -> Sing b -> Decision (a :<: b)
  infix 4 %<

-- | Instance for Peano Natural numbers. Uses unsafeCoerce function.
instance LSDecide Nat where
  (%<) a b | Just l <- lessNat a b = Proved l
           | otherwise = Disproved (\_ -> error "the first number is not less then the second")

lessNat :: forall (a :: Nat) (b :: Nat). Sing a -> Sing b -> Maybe (a :<: b)
lessNat a b | toNatural a <= toNatural b = Just (unsafeCoerce# (LS :: Zero :<: One))
            | otherwise = Nothing




-- | Without this lemma GHC can't deduce that (x <= x)
lemmaXLEqX :: forall (n :: Nat) (x :: Nat). (n ~ x) => Sing n -> Sing x -> (x :<: x)
lemmaXLEqX (SZ) (SZ) = LS
lemmaXLEqX (SS l) (SS k) = case lemmaXLEqX l k of LS -> LS

-- | Without this lemma GHC can' deduce that if ('S n <= x) then n <= x
lemma1 :: forall (n :: Nat) (x :: Nat). (('S n <= x) ~ 'True) => Sing n -> Sing x -> (n :<: x)
lemma1 SZ _ = LS
lemma1 (SS l) (SS k) = case lemma1 l k of LS -> LS
