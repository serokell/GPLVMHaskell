{-# LANGUAGE AllowAmbiguousTypes, EmptyCase, FlexibleContexts, LambdaCase,
             MagicHash, MultiParamTypeClasses, PolyKinds,
             UndecidableInstances #-}

module GPLVM.TypeSafe.PCA
       ( PCA (..)
       , covariance
       , eigenValues
       , eigenVectors
       , inputData
       , finalData
       , makePCA
       , restoredData
       , meanMatrix
       ) where

import Universum hiding (All, Nat, One, Vector, transpose)

import GPLVM.Types
import GPLVM.TypeSafe.Types
import GPLVM.TypeSafe.Util

import Control.Lens (makeLenses)
import Data.Array.Repa (D, extent)
import qualified Data.Array.Repa as R (Z)
import Data.Array.Repa.Repr.ForeignPtr
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

import Data.Singletons.Decide ((:~:)(..), Decision(..), (%~))
import Data.Type.Natural
import Data.Vinyl.TypeLevel (AllConstrained)

data PCA = PCA
  { _inputData    :: Matrix D Double
  , _covariance   :: Matrix D Double
  , _eigenVectors :: Matrix D Double
  , _eigenValues  :: Matrix D Double
  , _finalData    :: Matrix D Double
  , _restoredData :: Matrix D Double
  , _meanMatrix   :: Matrix D Double
  }

makeLenses ''PCA

data TypeSafePCA = forall x y d. (d <= x ~ 'True) => TypeSafePCA
  {
    desiredDim  :: Proxy d
  , inputData_  :: DimMatrix D y x Double
  , covariance_ :: DimMatrix D x x Double
  , eigenVectors_ :: DimMatrix D x x Double
  , eigenValues_ :: DimMatrix D x One Double
  , finalData_ :: DimMatrix D y d Double
  , restoredData_ :: DimMatrix D y x Double
  , meanMatrix_ :: DimMatrix D y x Double
  }


makePCA
    :: Int
    -> Matrix D Double
    -> PCA
makePCA desiredDimensions input =
  case toSing (intToNat desiredDimensions) of
    SomeSing (sd :: Sing desired) -> withSingI sd $ withMat input $ \(inputMatrix :: DimMatrix D y x Double) ->
      case checkInput (Proxy @desired) (Proxy @x) of
       Proved LS -> convertTypSafeToPCA $ makePCATypeSafe (Proxy @desired) inputMatrix
       Disproved _ -> error "t"

checkInput
  :: forall (d :: Nat) (x :: Nat). (SingI d, SingI x)
  => Proxy d
  -> Proxy x
  -> Decision (d :<: x)
checkInput _ _ =
  let des = (sing :: Sing d) %< (sing :: Sing x)
  in des

makePCATypeSafe
    :: forall d x y. (AllConstrained SingI '[d,x,y], d <= x ~ 'True)
    => Proxy (d :: Nat)
    -> DimMatrix D y x Double
    -> TypeSafePCA
makePCATypeSafe desiredDim inputData_ =
  let (meanVec, covariance_) = meanCovSM inputData_
      meanMatrix_ = extendY @y meanVec :: DimMatrix D y x Double
      adjustInput' = inputData_ -^^ meanMatrix_
      (eigenValues_, eigenVectors_) = eigSHM covariance_
      (eigenVecsf@(DimMatrix m)) = getFirstColumns @d eigenVectors_
      (finalData'@(DimMatrix m1)) = (transposeM eigenVecsf) `mulM` transposeM adjustInput' -- data in new eigenvectors space
      finalData_ = transposeM finalData'
      eigenVecsf' = transposeM eigenVecsf
      restoredDataWOMean@(DimMatrix m2) = trace (show @String $ extent m1) $ transposeM $ pinvSM eigenVecsf' `mulM` finalData'  -- restore to the original space
      restoredData_ = trace (show @String $ extent m2) $ restoredDataWOMean +^^ meanMatrix_              -- add mean
  in TypeSafePCA{..}

convertTypSafeToPCA :: TypeSafePCA -> PCA
convertTypSafeToPCA (TypeSafePCA _ input cov eigenVecs eigenVals final restored mean) =
  let (DimMatrix _inputData) = input
      (DimMatrix _covariance) = cov
      (DimMatrix _eigenVectors) = eigenVecs
      (DimMatrix _eigenValues) = eigenVals
      (DimMatrix _finalData) = final
      (DimMatrix _restoredData) = restored
      (DimMatrix _meanMatrix) = mean
  in PCA{..}
