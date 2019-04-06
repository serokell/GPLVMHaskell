module GPLVM.Util
       ( deleteRows
       , diag
       , getColumns
       , deleteColumns
       , getRows
       , mapDiagonal
       , mapIndexed
       , meanColumn
       , meanColumnWithNan
       , normalDistributionProbability
       , substractMeanWithNan
       , randomMatrixD
       , sumListMatrices
       , substractMean
       , toLatentSpacePoints
       , transposeMatrix
       , updateMatrix
       , updateMatrix2
       , toMatrix
       , toMatrix'
       , trace2S
       , trace2P
       , zipWithArray
       , (++^)
       , (**^)
       , (--^)
       , (//^)
       , (^++)
       , (^**)
       , (^--)
       , (^//)
       , (^--)
       , argmax
       ) where

import Prelude (RealFloat, head, isNaN)
import Universum hiding
  (All, Vector, head, map, size, sum, toList, trace, transpose, traverse,
  zipWith)

import Data.Array.Repa hiding ((++))
import Data.Array.Repa.Repr.Unboxed (Unbox)
import Data.Array.Repa.Shape (size)
import Data.Array.Repa.Unsafe (unsafeBackpermute)
import Data.List ((!!), (\\))
import Data.Random.Normal (normals)
import Debug.Trace
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector, diag)
import System.Random (Random, RandomGen)

import GPLVM.Types

trace2S :: Array U DIM2 Double -> Double
trace2S x
 = sumAllS $ unsafeBackpermute (Z :. (min nRows nColumns)) (\(Z :. i) -> (Z :. i :. i)) x
 where
    (Z :. nRows :. nColumns) = extent x

trace2P :: Monad m => Array U DIM2 Double -> m Double
trace2P x
 = sumAllP $ slice y (Z :. (0 :: Int) :. All)
 where
    y               =  unsafeBackpermute dim f x
    f (Z :. i :. j) = Z :. (i + j) `mod` nRows:. j
    dim@(Z :. nRows :. _) = extent x

deleteRows
  :: forall a. Unbox a
  => [Int]
  -> Matrix D a
  -> Matrix D a
deleteRows indices m@(ADelayed (Z :. y0 :. _) _) =
  getRows ([0..(y0 - 1)] \\ indices) m

getRows
  :: forall a. Unbox a
  => [Int]
  -> Matrix D a
  -> Matrix D a
getRows [] _ = delay $ fromListUnboxed (Z :. (0 :: Int) :. (0 :: Int)) []
getRows indices m@(ADelayed (Z :. y0 :. x0) _) =
  let numOfRows = length indices
  in if minimum indices < 0 || maximum indices > (y0 -1)
     then error $ "there are nonexistent numbers of rows"
     else backpermute (Z :. numOfRows :. x0) (\(Z :. y :. x) -> (Z :. (indices !! y) :. x)) m

getColumns
  :: forall a. Unbox a
  => [Int]
  -> Matrix D a
  -> Matrix D a
getColumns [] _ = delay $ fromListUnboxed (Z :. (0 :: Int) :. (0 :: Int)) []
getColumns indices m@(ADelayed (Z :. y0 :. x0) _) =
  let numOfColumns = length indices
  in if minimum indices < 0 || maximum indices > (x0 -1)
     then error "there are nonexistent numbers of columns"
     else backpermute (Z :. y0 :. numOfColumns) (\(Z :. y :. x) -> (Z :. y :. (indices !! x) )) m

deleteColumns
  :: forall a. Unbox a
  => [Int]
  -> Matrix D a
  -> Matrix D a
deleteColumns indices m@(ADelayed (Z :. _ :. x0) _) =
  getColumns ([0..(x0 - 1)] \\ indices) m

diag
  :: Matrix D a
  -> Matrix D a
diag m@(ADelayed (Z :. x0 :. y0) _) =
  if (x0 == y0)
  then fromFunction (Z :. x0 :. 1) (\(Z :. y :. _) -> m ! (Z :. y :. y))
  else error "should be squared !"

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

meanColumnWithNan
    :: forall a. (Unbox a, RealFloat a)
    => Matrix D a
    -> Vector D a
meanColumnWithNan matrix@(ADelayed (Z :. (n :: Int) :. _ ) _) =
  let sumVector = foldS (\acc x -> if (isNaN x) then acc else acc + x) 0.0 $ matrix
      lengthVector = foldS (\acc x -> if (isNaN x) then acc else acc + 1.0) 0.0 $ matrix
  in zipWith (/) sumVector lengthVector

substractMean :: forall a. (Unbox a, Fractional a)
    => Matrix D a
    -> Matrix D a
substractMean matrix = matrix ^-- (meanColumn matrix)

substractMeanWithNan :: forall a. (Unbox a, RealFloat a)
    => Matrix D a
    -> Matrix D a
substractMeanWithNan matrix = matrix ^-- (meanColumnWithNan matrix)


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

argmax
  :: (Ord a, Ord b)
  => (a -> b)
  -> [a]
  -> a
argmax f args =
  fst $
  maximum $
  [ (x, f x) | x <- args ]

normalDistributionProbability
  :: Double
  -> Double
  -> Double
  -> Double
normalDistributionProbability mean variance x =
  (1.0 / (sqrt $ 2*pi*variance)) * (exp $ ( (-1) * (x - mean) ^ 2) / (2*variance))
