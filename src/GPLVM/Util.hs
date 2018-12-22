module GPLVM.Util
       ( randomMatrixD
       , toLatentSpacePoints
       , transposeMatrix
       , updateMatrix
       , updateMatrix2
       ) where

import Universum hiding (transpose)

import Data.Array.Repa
import Data.Array.Repa.Repr.Unboxed (Unbox)

import Data.Random.Normal
import Numeric.LinearAlgebra.Repa hiding (Vector, Matrix)
import System.Random (Random, RandomGen)

import GPLVM.Types

updateMatrix
    :: (Array r DIM2 a -> Array r DIM2 a)
    -> Matrix r a
    -> Matrix r a
updateMatrix f (Matrix m) = Matrix (f m)

transposeMatrix
    :: Matrix D a
    -> Matrix D a
transposeMatrix = updateMatrix transpose

updateMatrix2
    :: (Array r DIM2 a -> Array r DIM2 a -> Array r DIM2 a)
    -> Matrix r a
    -> Matrix r a
    -> Matrix r a
updateMatrix2 f (Matrix m) (Matrix n) = Matrix (f m n)


toLatentSpacePoints
    :: forall a. Matrix D a
    -> LatentSpacePoints a
toLatentSpacePoints = LatentSpacePoints . updateMatrix transpose

randomMatrixD
    :: forall a g.
    ( RandomGen g
    , Random a
    , Unbox a
    , Floating a
    )
    => g
    -> (Int, Int)
    -> Matrix D a
randomMatrixD gen (rows, cols) =
    let randomList = take (rows * cols) (normals gen) in
    Matrix $ (smap id $ fromListUnboxed (Z :. rows :. cols) randomList)
