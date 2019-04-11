-- | based on Informative vector machine implemented as follows:
-- https://github.com/lawrennd/ivm/blob/master/matlab/ivmCreate.m
-- See also
-- https://papers.nips.cc/paper/2240-fast-sparse-gaussian-process-methods-the-informative-vector-machine.pdf

{-# LANGUAGE AllowAmbiguousTypes #-}
{-# LANGUAGE LambdaCase #-}

module GPLVM.TypeSafe.IVM
  ( IVM (..)
  , IVMTypeSafe (..)
  , Purpose (..)
  , makeIVM
  , makeIVMTypeSafe
  ) where

import Prelude (log, (!!))
import Universum hiding (All, Nat, One,
                         Vector, toList, transpose,
                         (%~))

import           Data.Array.Repa hiding (Z, zipWith)
import           Data.Array.Repa (toList)
import           Data.Default (Default (def))
import           Data.Singletons.Decide (Decision (..), (:~:) (..), (%~))
import           Data.Type.Natural
import           Data.Vinyl.TypeLevel (AllConstrained)

import           GPLVM.Types (Matrix (..))
import           GPLVM.TypeSafe.Types (DimMatrix (..), (:<:) (..), (%<))
import           GPLVM.TypeSafe.Util (withMat, Decisions)
import           GPLVM.Util (argmax, normalDistributionProbability)

import           Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

import           Statistics.Distribution (cumulative)
import           Statistics.Distribution.Normal (standard)

data IVMTypeSafe =
  forall (d :: Nat)                            -- ^ d is a number of active points
  (m :: Nat)                                   -- ^ input columns
  (n :: Nat).                                  -- ^ input rows
  (d <= n ~ 'True) => IVMTypeSafe              -- ^ a justification that d <= n is true
  { covariance    :: DimMatrix D n n Double    -- ^ input covariance
  , inputMatrix   :: DimMatrix D m n Double    -- ^ input matrix
  , sparsedMatrix :: DimMatrix D m d Double    -- ^ output sparsed matrix
  , resultIndices :: [Int]                     -- ^ the list of row indices in original
                                               -- | that were placed to the output sparsed matrix
  , bias          :: Double                    -- ^ bias parameter
  }

-- dimesionless version of IVM
data IVM = IVM
  { covar         :: Matrix D Double  -- ^ dimensionless covariance
  , input         :: Matrix D Double  -- ^ dimesionless input matrix
  , sparsed       :: Matrix D Double  -- ^ dimesionless sparsed output matrix
  , selectedRows  :: [Int]            -- ^ selected rows defined similarly to resultIndices
                                      -- | in IVMTypeSafe
  , bias'         :: Double           -- ^ bias is like bias in IVMTypeSafe
  }

-- | We may use IVM either for regression or classification.
-- | This flag allows one to choose one of these tasks
data Purpose =
    Regression
  | Classification
  deriving (Show, Read)


-- | By default, IVM makes regression
instance Default Purpose where
  def = Regression

makeIVMTypeSafe
  :: forall (d :: Nat) (m1 :: Nat) (m2 :: Nat) (n1 :: Nat) (n2 :: Nat).
  ( AllConstrained SingI [m1, m2, n1, n2, d]
  , d <= n1 ~ 'True
  , m1 ~ n1
  , n1 ~ n2
  )
  => DimMatrix D m1 n1 Double  -- ^ covariance
  -> DimMatrix D m2 n2 Double  -- ^ input
  -> Double                    -- ^ bias parameter
  -> Purpose                   -- ^ purpose: either regression or classification
  -> IVMTypeSafe
makeIVMTypeSafe cov input biasP = \case
  Regression -> actionReg
  Classification -> actionClass
  where
    activeNum, actualNum
      :: Int
    activeNum = natToInt $ demote @d :: Int
    actualNum = natToInt $ demote @m2 :: Int

    actionReg  = undefined

    actionClass = undefined

    allIndices = [1..actualNum]

    activeIndices = [1..activeNum]

    internalStep xs ys =
      let l1 = zipWith gIN xs ys
          l2 = zipWith nuIN xs ys
          l3 = zipWith deltaH xs ys
      in (l1, l2, l3)

    yN :: Int -> Double
    yN = undefined

    -- cumulative Gaussian distribution
    phi :: Double -> Double
    phi = cumulative standard

    -- | the rest formulae required for IVM iterations
    gIN, nuIN, muIN, uIN, ksiIN, deltaH
      :: Int -> Int -> Double

    muIN i n = undefined

    uIN i n =
      cIN i n * muIN i n + biasP

    gIN i n = (cIN (i - 1) n) +
      (normalDistributionProbability 0 1 (uIN (i - 1) n)) +
    	(phi $ uIN (i - 1) n)

    cIN i n = yN n * (sqrt (ksiIN i n))

    nuIN i n =
      undefined

    ksiIN i n =
      take i (toList covarianceDiag) !! n

    -- | vector obtained from covariance diagonal
    covarianceDiag =
      diagD . takeDiagD $ getInternal cov

    -- | differential entropy
    deltaH i n =
      - 0.5 * (log $ 1 - nuIN i n * ksiIN (i - 1) n)

    -- vectors that should be updated during the main iteration

    argMaxDelta i xs = argmax $ deltaH i <$> xs

    sIN i n = undefined



makeIVM
  :: Int              -- ^ a number of active point
  -> Matrix D Double  -- ^ covariance matrix
  -> Matrix D Double  -- ^ input matrix
  -> Double           -- ^ bias parameter
  -> Purpose
  -> IVM
makeIVM actPoints cov input biasP purpose =
  case toSing (intToNat actPoints) of
    SomeSing (active :: Sing active) -> withSingI active $ withMat cov $
      \(mat1 :: DimMatrix D y x Double) -> withMat input $
      \(mat2 :: DimMatrix D y2 x2 Double) ->
      let (sol1, sol2, sol3) = checkInput @y @x @y2 @x2 @active in
      case (sol1, sol2) of
        _ -> error "equalities are false"
        (Proved Refl, Proved Refl) -> case sol3 of
          Disproved _ -> error "desired dimension is greater than required"
          Proved LS   -> convertTypeSafeIVM $ makeIVMTypeSafe @active mat1 mat2 biasP purpose
  where
    checkInput
      :: forall (y :: Nat) (x :: Nat) (y2 :: Nat) (x2 :: Nat) (active :: Nat).
      (AllConstrained SingI '[x, y, x2, active])
      => Decisions (x :~: y, x :~: x2, active :<: x)
    checkInput =
      let des1 = (sing :: Sing x) %~ (sing :: Sing y)
          des2 = (sing :: Sing x) %~ (sing :: Sing x2)
          des3 = (sing :: Sing active) %< (sing :: Sing x)
      in (des1, des2, des3)

    convertTypeSafeIVM
      :: IVMTypeSafe
      -> IVM
    convertTypeSafeIVM (IVMTypeSafe (DimMatrix cov) (DimMatrix input) (DimMatrix sparsed) ind bias) =
      IVM cov input sparsed ind bias
