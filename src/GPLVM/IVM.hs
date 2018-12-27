{-
based on Informative vector machine implemented based on the following:
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.440.3620&rep=rep1&type=pdf

-}

module GPLVM.IVM where

import Universum hiding (Vector)

import Control.Lens (makeLenses)
import Data.Array.Repa
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

import GPLVM.Types

newtype IVMInput a = IVMInput { unInput :: Matrix a }

newtype IVMOutput a = IVMOutput { unOutput :: Matrix a }

----- simple ivm

simpleIVM
    :: IVMInput a
    -> Int      --- a new number of rows
    -> IVMOutput a
simpleIVM input@(IVMInput matrixInput) =



----- multitask ivm
