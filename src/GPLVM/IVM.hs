-- | based on Informative vector machine implemented as follows:
-- https://github.com/lawrennd/ivm/blob/master/matlab/ivmCreate.m
-- See also
-- https://papers.nips.cc/paper/2240-fast-sparse-gaussian-process-methods-the-informative-vector-machine.pdf


module GPLVM.IVM
       ( IVMInput (..)
       , IVMOutput (..)
       , runIVM
       ) where

import Universum hiding (Vector)

import Prelude (log)

import Control.Lens (makeLenses)
import Data.Array.Repa
import GPLVM.TypeSafe.Types
import GPLVM.TypeSafe.Util
--import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
