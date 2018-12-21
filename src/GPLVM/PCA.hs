module GPLVM.PCA where

import Universum hiding (Vector)

import GPLVM.Types

import Control.Monad.Identity
import Data.Array.Repa hiding (Vector)
import Data.Array.Repa.Repr.ForeignPtr
import Data.Bifunctor
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import Numeric.LinearAlgebra.Repa.Conversion
import qualified Numeric.Matrix as M

data PCA = PCA
  {   inputData :: Matrix D Double
    , covariance  :: Covariance Double
    , eigenVectors :: EigenVectors Double
    , eigenValues :: Vector F Double
    , finalData :: Matrix D Dbouble
    , restoredData :: Matrix D Double 
  }

makePCA :: Int -> Matrix D Double -> PCA
makePCA desiredDimentions input =
  let inputData@(Matrix inp) = input
      (Z :. yInp :. xInp) = extent inp 
      (meanVec, covariance) = runIdentity $ meanCovP inp
      adjustInput = inp -^ (extend (Z :. yInp :. all) meanVec)
      (eigenValues, eigenVec) = first Vector $ eigSH covariance
      (AForeignPtr (Z :. y :. x) _ _) = eigenVec
      eigenVectors =
        if desiredDimentions >= y
        then Matrix eigen
        else Matrix $ fromForeignPtr (Z :. desiredDimentions :. x) $ toForeignPtr eigenVec
      finalData = eigenVectors `mul` adjustInput
      restoredData = inv eigenVectors `mul` finalData
  in PCA{..}

