{-
based on Informative vector machine implemented as follows:
https://github.com/lawrennd/ivm/blob/master/matlab/ivmCreate.m
-}

module GPLVM.IVM where

import Universum hiding (Vector)

import Control.Lens (makeLenses)
import Data.Array.Repa
import GPLVM.TypeSafe.Types
import GPLVM.TypeSafe.Util
--import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

