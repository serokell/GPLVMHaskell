module GPLVM.Util
       ( mapDiagonal
       , mapIndexed
       , meanColumn
       , randomMatrixD
       , sumListMatrices
       , substractMean
       , toLatentSpacePoints
       , transposeMatrix
       , updateMatrix
       , updateMatrix2
       , (^--)
       ) where

import Prelude (head)
import Universum hiding (Vector, head, map, size, transpose, traverse, zipWith)

import Data.Array.Repa
import Data.Array.Repa.Repr.Unboxed (Unbox)
import Data.Array.Repa.Shape (size)

import Data.Random.Normal
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import System.Random (Random, RandomGen)

import GPLVM.Types

updateMatrix
    :: (Array r DIM2 a -> Array r DIM2 a)
    -> Matrix r a
    -> Matrix r a
updateMatrix f m = f m

transposeMatrix
    :: Matrix D a
    -> Matrix D a
transposeMatrix = updateMatrix transpose

updateMatrix2
    :: (Array r DIM2 a -> Array r DIM2 a -> Array r DIM2 a)
    -> Matrix r a
    -> Matrix r a
    -> Matrix r a
updateMatrix2 f m n = f m n


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
    delay $ fromListUnboxed (Z :. rows :. cols) randomList

sumListMatrices
    :: forall a. (Unbox a, Num a)
    => [Matrix D a]
    -> Matrix D a
sumListMatrices = foldl1 (\acc m -> acc +^ m)

meanColumn
    :: forall a. (Unbox a, Fractional a)
    => Matrix D a
    -> Vector D a
meanColumn matrix@(ADelayed (Z :. (n :: Int) :. _ ) _) =
  map (/(fromIntegral n)) . foldS (+) 0 $ transpose matrix

substractMean :: forall a. (Unbox a, Fractional a)
    => Matrix D a
    -> Matrix D a
substractMean matrix = matrix ^-- (meanColumn matrix)

mapIndexed
    :: forall a b sh r. (Unbox a, Unbox b, Shape sh, Source r a)
    => (a -> sh -> b)
    -> Array r sh a
    -> Array D sh b
mapIndexed f array = traverse array id (\g sh -> f (g sh) sh)

mapDiagonal
  :: forall a sh r. (Unbox a, Shape sh, Source r a)
  => (a -> a)
  -> Array r sh a
  -> Array D sh a
mapDiagonal f = mapIndexed (\value sh -> if (listEquality $ listOfShape sh) then f value else value)

listEquality
    :: forall a. Eq a
    => [a]
    -> Bool
listEquality [] = True
listEquality lst =
  let first = head lst
  in all (\x -> x == first) lst

toMatrix
    :: Source r1 a
    => Array r1 DIM1 a
    -> Int
    -> Array D DIM2 a
toMatrix arr desiredSize =
    fromFunction (Z :. desiredSize :. newDimension) generator
    where
        generator (Z :. rows :. cols) = linearIndex arr cols
        newDimension = size $ extent $ arr

toMatrix'
    :: Source r1 a
    => Array r1 DIM1 a
    -> Array D DIM2 a
toMatrix' arr =
    fromFunction (Z :. dimension :. 1) generator
    where
        dimension = size . extent $ arr
        generator (Z :. rows :. cols) = linearIndex arr rows

zipWithArray
    :: (a -> b -> c)
    -> Array D DIM1 a
    -> Array D DIM2 b
    -> Array D DIM2 c
zipWithArray f array1 array2@(ADelayed (Z :. (n :: Int) :. _) _) =
    zipWith f (toMatrix array1 n) array2

zipWithArray'
    :: (a -> b -> c)
    -> Array D DIM2 a
    -> Array D DIM1 b
    -> Array D DIM2 c
zipWithArray' f array1@(ADelayed (Z :. (n :: Int) :. _) _) array2 =
    zipWith f array1 (toMatrix array2 n)

infixl 6 ++^
infixl 6 --^
infixl 7 **^
infixl 7 //^

(++^), (**^), (--^), (//^)
    :: (Num a, Fractional a)
    => Array D DIM1 a
    -> Array D DIM2 a
    -> Array D DIM2 a
(++^) = zipWithArray (+)
(**^) = zipWithArray (*)
(--^) = zipWithArray (-)
(//^) = zipWithArray (/)

infixl 6 ^++
infixl 6 ^--
infixl 7 ^**
infixl 7 ^//


(^++), (^**), (^--), (^//)
    :: (Num a, Fractional a)
    => Array D DIM2 a
    -> Array D DIM1 a
    -> Array D DIM2 a
(^++) = zipWithArray' (+)
(^**) = zipWithArray' (*)
(^--) = zipWithArray' (-)
(^//) = zipWithArray' (/)
