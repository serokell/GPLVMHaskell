module GPLVM.Basic where

import Universum hiding (transpose)

import           Control.Lens (makeLenses)

import           Data.Vector (Vector)
import           Data.Vector as V
-- import qualified Data.Sequence as S
import           Numeric.Matrix (Matrix, MatrixElement)
import qualified Numeric.Matrix as M

import           GPLVM.LinearSolver
import           GPLVM.Util

newtype EigenValue a = EigenValue { _unEigen :: a }

makeLenses ''EigenValue

data EigenValueException =
      GivenValueNotEigen
    | GivenVectorNotEigen
    | OutOfLength
    | LinearDependency

toEigenValue
     :: forall a. MatrixElement a
     => Matrix a
     -> EigenValue a
     -> Either EigenValueException (EigenValue a)
toEigenValue matrix v@(EigenValue value) = do
     let matrixRows = M.numRows matrix
     let nUnit = M.unit matrixRows
     let eigenDet = M.det (matrix `M.minus` (nUnit `M.scale` value))
     case (M.isSquare matrix && eigenDet == 0) of
         True -> return v
         False -> Left GivenValueNotEigen

isEigenVector
    :: forall a. MatrixElement a
    => Matrix a
    -> EigenValue a
    -> Vector a
    -> Either EigenValueException (Vector a)
isEigenVector matrix v vec = do
    val <- toEigenValue matrix v
    let vecMat = vectorToMatrix vec
    let firstPredicate = V.length vec == M.numCols matrix
    let secondPredicate = matrix `M.times` vecMat == vecMat `M.scale` (val ^. unEigen)
    case firstPredicate && secondPredicate of
        False -> Left OutOfLength
        True -> return vec


toEigenVector ::
    forall a. MatrixElement a
    => Matrix a
    -> EigenValue a
    -> LinearSolver a
    -> Either EigenValueException (Vector a)
toEigenVector matrix v solver = do
    val <- toEigenValue matrix v
    let value = val ^. unEigen
    case M.det matrix == 0 of
        False -> Left LinearDependency
        True -> do
            let matrixRows = M.numRows matrix
            let nUnit = M.unit matrixRows
            let zeroColumn = V.fromList $ M.trace nUnit
            let matrix' = matrix `M.minus` (nUnit `M.scale` value)
            return $ solver matrix' zeroColumn

index
    :: forall a. MatrixElement a
    => (Int, Int)
    -> Matrix a
    -> Maybe a
index (rows, cols) matrix = do
    case (rows < M.numRows matrix, cols < M.numCols matrix) of
        (True, True) -> return $ M.at matrix (rows, cols)
        _ -> Nothing

{-
replace
    :: forall a. MatrixElement a
    => (Int, Int)
    -> a
    -> Matrix a
    -> Maybe (Matrix a)
replace (givenRow, givenCol) element matrix = undefined

choleskyDecomposition
    :: forall a. MatrixElement a
    => Matrix a
    -> Matrix a
choleskyDecomposition matrix = undefined -}
