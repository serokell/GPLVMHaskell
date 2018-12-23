module GPLVM.GaussianProcessLV
       ( GaussianProcessLatentVariable (..)
       , gP
       , initializeLatent
       , lvObserved
       ) where

import Universum

import Control.Lens (makeLenses)

import Data.Array.Repa

import GPLVM.GaussianProcess
import GPLVM.PCA
import GPLVM.Types (Matrix (..), ObservedData (..), unObservedData)

data GaussianProcessLatentVariable = GaussianProcessLatentVariable
    { _lvObserved :: ObservedData Double
    , _gP :: GaussianProcess Double
    }

makeLenses ''GaussianProcessLatentVariable

initializeLatent
    :: Int
    -> ObservedData Double
    -> Matrix D Double
initializeLatent latentDimension observedData =
    (makePCA latentDimension (observedData ^. unObservedData)) ^. finalData


{-
toGPLV
    :: Int
       -- latent dimension
    -> ObservedData Double
       -- observed data
    -> GaussianProcessLatentVariable
toGPLV latentDimension observedData =
    initializeLatent latentDimension observedData
-}
