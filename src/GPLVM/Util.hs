module GPLVM.Util where

import           Universum

import           Data.Vector (Vector)
import qualified Data.Vector as V

import           Numeric.Matrix (Matrix, MatrixElement)
import qualified Numeric.Matrix as M

foldrMatrix
    :: MatrixElement e =>
    (e -> b -> b) -> b -> Matrix e -> b
foldrMatrix f element matrix =
    appEndo (M.foldMap (Endo . f) matrix) element

vectorToMatrix
    :: forall a. MatrixElement a
    => Vector a -> Matrix a
vectorToMatrix vec = M.transpose . M.fromList $ [V.toList vec]

vectorListToMatrix
    :: forall a. MatrixElement a
    => [Vector a] -> Matrix a
vectorListToMatrix vecs = M.fromList $ V.toList <$> vecs
