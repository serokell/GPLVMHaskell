module Math.TestDimensionsPCA where

import Universum hiding (All, Vector, transpose)

import GPLVM.Types

import Control.Lens (makeLenses)
import Data.Array.Repa
import Data.Array.Repa.Repr.ForeignPtr
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

import Math.Matrix

data PCA = PCA
  { inputData    :: Matrix D Double
  , covariance   :: Covariance Double
  , eigenVectors :: EigenVectors Double
  , eigenValues  :: Vector F Double
  , finalData    :: Matrix D Double
  , restoredData :: Matrix D Double
  }

makePCA
    :: Int
    -> Matrix D Double
    -> PCA
makePCA desiredDimensions input =
  let inputData = input
      (Z :. yInp :. _) = extent input                                          -- get dimension
      (meanVec, covariance) = meanCovS input
      meanMatrix = DimMatrix $ toMatrix meanVec yInp
      adjustInput = DimMatrix input -^^ meanMatrix
      (eigenValues', eigenVec) = eigSHM covariance
      eigenValue = runVector eigenValues'
      eigenVecD@(DimMatrix (ADelayed (Z :. y :. x) f)) = transposeM eigenVec              -- transposing is almost free
      eigenVectors' =                                                        -- leave only n desired eigenvectors
        if desiredDimensions >= y
        then eigenVecD
        else DimMatrix $ ADelayed (Z :. desiredDimensions :. x) f
      eigenVectors = getInternal $ transposeM eigenVectors'                       -- colunmns are eigenvectors
      finalData' = eigenVectors' `mulPM` transposeM adjustInput -- data in new eigenvectors space
      finalData = getInternal . transposeM $ finalData'
      restoredDataWOMean = transposeM $ pinvSM eigenVectors' `mulPM` finalData'  -- restore to the original space
      restoredData = getInternal $ restoredDataWOMean +^^ meanMatrix              -- add mean
  in PCA{..}
