module GPLVM.GPTest where

import Universum
import Prelude ((!!))

import Data.Array.Repa
import Numeric.LinearAlgebra.Repa hiding (Matrix)

import GPLVM.Types
import GPLVM.Util

testList :: [Double]
testList = [ -5.0
           , -4.79591837
           , -4.59183673
           , -4.3877551
           , -4.18367347
           , -3.97959184
           , -3.7755102
           , -3.57142857
           , -3.36734694
           , -3.16326531
           , -2.95918367
           , -2.75510204
           , -2.55102041
           , -2.34693878
           , -2.14285714
           , -1.93877551
           , -1.73469388
           , -1.53061224
           , -1.32653061
           , -1.12244898
           , -0.91836735
           , -0.71428571
           , -0.51020408
           , -0.30612245
           , -0.10204082
           ,  0.10204082
           ,  0.30612245
           ,  0.51020408
           ,  0.71428571
           ,  0.91836735
           ,  1.12244898
           ,  1.32653061
           ,  1.53061224
           ,  1.73469388
           ,  1.93877551
           ,  2.14285714
           ,  2.34693878
           ,  2.55102041
           ,  2.75510204
           ,  2.95918367
           ,  3.16326531
           ,  3.36734694
           ,  3.57142857
           ,  3.7755102
           ,  3.97959184
           ,  4.18367347
           ,  4.3877551
           ,  4.59183673
           ,  4.79591837
           ,  5.0
           ]

xTestDim :: DIM2
xTestDim = (Z :. 50 :. 1)

xTestFunction :: DIM2 -> Double
xTestFunction (Z :. n :. _) = testList !! n

xTest :: Matrix D Double
xTest = fromFunction xTestDim xTestFunction

inputObserve :: InputObservations Double
inputObserve = InputObservations xTest

foo :: Array D DIM2 Double -> Array U DIM2 Double
foo = computeS @D @_ @_ @U

matrixSquare :: Matrix D Double -> Matrix D Double
matrixSquare matrix = smap ((flip (^)) 2) matrix

fun :: Double -> Double
fun = exp . ((*) (-0.5)) . ((*) 10)

sqDist :: Matrix D Double -> Matrix D Double -> Matrix D Double
sqDist matrix1 matrix2 = matrixSquare matrix1 +^ matrix2'
     where
       matrix2' = matrixSquare matrix2 -^
                  (delay $ smap ((*) 2) $ matrix1 `mulS` (transposeMatrix matrix2))

kernelFunction :: KernelFunction Double
kernelFunction matrix1 matrix2 =
    smap fun $ sqDist matrix1 matrix2

xTrainList :: [Double]
xTrainList = [-4, -3, -2, -1, 1]

xTrainFunction :: DIM2 -> Double
xTrainFunction (Z :. n :. _) = xTrainList !! n

xTrain :: Matrix D Double
xTrain = fromFunction (Z :. 5 :. 1) xTrainFunction

yTrain :: Matrix D Double
yTrain = smap sin xTrain
