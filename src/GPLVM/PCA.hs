module GPLVM.PCA
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

import Universum hiding (All, Vector, transpose)

import GPLVM.Types

import Control.Lens (makeLenses)
import Data.Array.Repa
import Data.Array.Repa.Repr.ForeignPtr
import Debug.Trace
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

data PCA = PCA
  {   _inputData    :: Matrix D Double
    , _covariance   :: Covariance Double
    , _eigenVectors :: EigenVectors Double
    , _eigenValues  :: Vector F Double
    , _finalData    :: Matrix D Double
    , _restoredData :: Matrix D Double
    , _meanMatrix   :: Matrix D Double
  }

makeLenses ''PCA

makePCA
    :: Int
    -> Matrix D Double
    -> PCA
makePCA desiredDimensions input =
  let _inputData = input
         -- get dimension
      (Z :. yInp :. _) = extent input                                          
      (meanVec, _covariance) = meanCovS input
      _meanMatrix = extend (Z :. yInp :. All) meanVec
      adjustInput'@(ADelayed (Z :. y1 :. x1) _) = input -^ _meanMatrix
      (_eigenValues, eigenVec) = eigSH _covariance
      eigenVecD@(ADelayed (Z :. y :. x) f) = transpose eigenVec        
        -- leave only n desired eigenvectors
      eigenVectors' =                                                        
        if desiredDimensions >= y
        then eigenVecD
        else ADelayed (Z :. desiredDimensions :. x) f
        -- colunmns are eigenvectors
      _eigenVectors = transpose eigenVectors'                      
        -- data in new eigenvectors space
      finalData' = runIdentity $ eigenVectors' `mulP` transpose adjustInput' 
      _finalData = transpose finalData'
        -- restore to the original space
      restoredDataWOMean = transpose $ pinvS eigenVectors' `mul` finalData'
        -- add mean  
      _restoredData = restoredDataWOMean +^ _meanMatrix              
  in PCA{..}
