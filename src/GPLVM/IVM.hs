{-
based on Informative vector machine implemented based on the following:
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.440.3620&rep=rep1&type=pdf

-}

module GPLVM.IVM
       ( IVMInput (..)
       , IVMOutput (..)
       , simpleIVM
       , unInput
       , unOutput
       ) where

import Universum hiding (Vector)

import Control.Lens (makeLenses)
import Data.Array.Repa
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

import GPLVM.Types

newtype IVMInput a = IVMInput { unInput :: Matrix D a }

newtype IVMOutput a = IVMOutput { unOutput :: Matrix D a }

makeLenses ''IVMInput
makeLenses ''IVMOutput

----- simple ivm

simpleIVM
    :: IVMInput a
    -> Int      --- a new number of rows
    -> IVMOutput a
simpleIVM input@(IVMInput matrixInput) rows =
    undefined

----- multitask ivm
