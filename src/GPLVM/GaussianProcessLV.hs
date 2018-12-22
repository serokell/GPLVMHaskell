module GPLVM.GaussianProcessLV where

import Universum

import Control.Lens (makeLenses)

import Data.Array.Repa
import GPLVM.Types (Distribution, Matrix (..), ObservedData (..))

data GaussianProcessLatentVariable a = GaussianProcessLatentVariable {
      _GPLVMObserved :: ObservedData a
    , _kerGPLVM :: Matrix D a -> Matrix D a -> Matrix D a
    , _distrGPLVM ::  Distribution D a
    }

makeLenses ''GaussianProcessLatentVariable
