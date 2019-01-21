module GPLVM.GPExample where

import Universum hiding (transpose, Vector)

import Data.Vector.Unboxed as V (fromList)

import Data.Array.Repa
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

import GPLVM.GaussianProcess
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

xTest :: Vector D Double
xTest = delay $ fromUnboxed (Z :. 50) (V.fromList testList)

inputObserve :: InputObservations Double
inputObserve = InputObservations xTest

computeSD
    :: Shape sh
    => Array D sh Double
    -> Array U sh Double
computeSD = computeS @D @_ @_ @U

sqDist
    :: Vector D Double
    -> Vector D Double
    -> Matrix D Double
sqDist vec1 vec2 = matrixSquare vec1 ++^ matrix2'
     where
         matrix2' = matrixSquare vec2 --^
                      ((delay . smap ((*) 2)) $
                      matrix1 `mulS` matrix2)
         matrixSquare matrix = smap flipExp matrix
         flipExp = (flip (^)) 2
         matrix1 = toMatrix' $ vec1
         matrix2 = transposeMatrix . toMatrix' $ vec2

kernelFunction
    :: Vector D Double
    -> Vector D Double
    -> Matrix D Double
kernelFunction matrix1 matrix2 =
    smap fun $ sqDist matrix1 matrix2
    where
        fun = exp . ((*) (-0.5)) . ((*) (1/0.1))

xTrainList :: [Double]
xTrainList =
    [ -4
    , -3
    , -2
    , -1
    , 1
    ]

xTrain :: Vector D Double
xTrain = delay $ fromUnboxed (Z :. 5) (V.fromList xTrainList)

yTrain :: Vector D Double
yTrain = smap sin xTrain

testTrainingData :: GPTrainingData Double
testTrainingData = GPTrainingData xTrain yTrain

meanTest :: GPMean Double
meanTest matrix = undefined

testGaussianProcess :: GaussianProcess Double
testGaussianProcess = GaussianProcess kernelFunction meanTest
