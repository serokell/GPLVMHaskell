module GPLVM.LinearSolver
       ( LinearSolver
       , LinearSolverException
       ) where

import Universum

import           Data.List ((!!))
import qualified Data.List.NonEmpty as N
import           Data.Vector (Vector)
import qualified Data.Vector as V

import           Numeric.Matrix (Matrix, MatrixElement)
import qualified Numeric.Matrix as M

import           GPLVM.Util

data LinearSolverException = HasNoSolution

type LinearSolver a = Matrix a -> Vector a -> Vector a
