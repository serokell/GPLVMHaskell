module GPLVM.PCA where

import Universum hiding (Vector, All)

import GPLVM.Types

import Control.Lens (makeLenses)
import Data.Array.Repa
import Data.Array.Repa.Repr.ForeignPtr
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

data PCA = PCA
  {   _inputData :: Matrix D Double
    , _covariance  :: Covariance Double
    , _eigenVectors :: EigenVectors Double
    , _eigenValues :: Vector F Double
    , _finalData :: Matrix D Double
    , _restoredData :: Matrix D Double
  }

makeLenses ''PCA

makePCA
    :: Int
    -> Matrix D Double
    -> PCA
makePCA desiredDimensions input =
  let _inputData@(Matrix inp) = input
      (Z :. yInp :. _) = extent inp
      (meanVec, _covariance) = runIdentity $ meanCovP inp
      adjustInput = inp -^ (extend (Z :. yInp :. All) meanVec)
      (_eigenValues, eigenVec) = first Vector $ eigSH _covariance
      (AForeignPtr (Z :. y :. x) _ _) = eigenVec
      eigenVectors' =
        if desiredDimensions >= y
        then eigenVec
        else fromForeignPtr (Z :. desiredDimensions :. x) $ toForeignPtr eigenVec
      _eigenVectors = Matrix eigenVectors'
      finalData' =
          runIdentity $ smap id eigenVectors' `mulP` adjustInput
      _finalData = Matrix $ smap id $ finalData'
      _restoredData = Matrix $ smap id $ inv eigenVectors' `mul` finalData'
  in PCA{..}
