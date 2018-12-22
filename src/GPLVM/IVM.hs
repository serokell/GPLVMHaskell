{-
based on PPCA implemented as follows:
https://github.com/lawrennd/ivm/blob/master/matlab/ivmCreate.m
-}

module GPLVM.IVM
       ( PPCA (..)
       ) where

import Universum hiding (Vector)

import Control.Lens (makeLenses)
import Data.Array.Repa
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

import GPLVM.Types (Matrix (..), ObservedData (..), Vector (..))

data ResultsPPCA a = ResultsPPCA
    { _coeffMatrix :: Matrix D a           --W
    , _condEst :: Matrix D a               -- conditional expectation
    , _lookLikelihood :: a
    }

data PPCA a = PPCA
    { _observedPPCA :: ObservedData a
    , _coeffPPCA :: Matrix D a            -- principal components coefficient
    , _scorePPCA :: Matrix D a            -- principal component scores
    , _variancePC :: Matrix D a           -- principal component variances
    , _meanPPCA :: Vector D a
    , _outputPPCA :: ResultsPPCA a
    }

makeLenses ''ResultsPPCA
makeLenses ''PPCA

{-
makePPCA
    :: Int
    -> ObservedData Double
    -> PPCA
-}
