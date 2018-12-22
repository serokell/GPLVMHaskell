module GPLVM.PCA
       ( PCA (..)
       , makePCA
       ) where

import Universum hiding (All, Vector, transpose)

import GPLVM.Types

import Control.Lens (makeLenses)
import Data.Array.Repa
import Data.Array.Repa.Repr.ForeignPtr
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

data PCA = PCA
  {   _inputData    :: Matrix D Double
    , _covariance   :: Covariance Double
    , _eigenVectors :: EigenVectors Double
    , _eigenValues  :: Vector F Double
    , _finalData    :: Matrix D Double
    , _restoredData :: Matrix D Double
  }

makeLenses ''PCA

makePCA
    :: Int
    -> Matrix D Double
    -> PCA
makePCA desiredDimensions input =
  let _inputData@(Matrix inp) = input
      (Z :. yInp :. _) = extent inp                                          -- get dimension
      (meanVec, _covariance) = meanCovS inp
      meanMatrix = extend (Z :. yInp :. All) meanVec
      adjustInput' = inp -^ meanMatrix
      (_eigenValues, eigenVec) = first Vector $ eigSH _covariance
      eigenVecD@(ADelayed (Z :. y :. x) f) = transpose eigenVec              -- transposing is almost free
      eigenVectors' =                                                        -- leave only n desired eigenvectors
        if desiredDimensions >= y
        then eigenVecD
        else ADelayed (Z :. desiredDimensions :. x) f
      _eigenVectors = Matrix $ transpose eigenVectors'                       -- colunmns are eigenvectors
      finalData' = runIdentity $ eigenVectors' `mulP` transpose adjustInput' -- data in new eigenvectors space
      _finalData = Matrix $ transpose finalData'
      restoredDataWOMean = transpose $ pinvS eigenVectors' `mul` finalData'  -- restore to the original space
      _restoredData = Matrix $ restoredDataWOMean +^ meanMatrix              -- add mean
  in PCA{..}
