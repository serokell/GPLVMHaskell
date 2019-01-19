module GPLVM.GPTest where

import Universum hiding (transpose, Vector)
import Prelude ((!!))

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

xTest :: Matrix D Double
xTest = delay $ fromUnboxed (Z :. 50 :. 1) (V.fromList testList)

inputObserve :: InputObservations Double
inputObserve = InputObservations xTest

foo :: Array D DIM2 Double -> Array U DIM2 Double
foo = computeS @D @_ @_ @U

matrixSub
    :: Num a
    => Vector D a
    -> Matrix D a
    -> Matrix D a
matrixSub vec@(ADelayed (Z :. len) f) matrix@(ADelayed (Z :. rows :. cols) g) =
    case (rows == len) of
        False -> error "matrixSub: length of vector should be equal to the number of rows"
        True -> delay $ newMatrix -^ matrix
    where
        newMatrix = fromFunction (Z :. rows :. cols) subFunction
        subFunction (Z :. rows' :. cols') = f (Z :. rows')

matrixSum
    :: Num a
    => Vector D a
    -> Matrix D a
    -> Matrix D a
matrixSum vec@(ADelayed (Z :. len) f) matrix@(ADelayed (Z :. rows :. cols) g) =
    case (rows == len) of
        False -> error "matrixSub: length of vector should be equal to the number of rows"
        True -> delay $ newMatrix +^ matrix
    where
        newMatrix = fromFunction (Z :. rows :. cols) sumFunction
        sumFunction (Z :. rows' :. cols') = f (Z :. rows')

sqDist :: Matrix D Double -> Matrix D Double -> Matrix D Double
sqDist matrix1 matrix2 = (toVector . matrixSquare $ matrix1) `matrixSum` matrix2'
     where
         matrix2' = (toVector . matrixSquare $ matrix2) `matrixSub`
                      ((delay . smap ((*) 2)) $ matrix1 `mulS` (transposeMatrix matrix2))
         matrixSquare matrix = smap flipExp matrix
         flipExp = (flip (^)) 2

toVector :: Matrix D a -> Vector D a
toVector matrix@(ADelayed (Z :. rows :. cols) f) =
    case cols == 1 of
        True -> fromFunction (Z :. rows) (\(Z :. x) -> f (Z :. x :. 1))
        False -> error "col should be single"

kernelFunction :: KernelFunction Double
kernelFunction matrix1 matrix2 =
    smap fun $ sqDist matrix1 matrix2
    where
        fun = exp . ((*) (-0.5)) . ((*) (1/0.1))

xTrainList :: [Double]
xTrainList =
    [ -4.0
    ,-3.83673469
    ,-3.67346939
    ,-3.51020408
    ,-3.34693878
    ,-3.18367347
    ,-3.02040816
    ,-2.85714286
    ,-2.69387755
    ,-2.53061224
    ,-2.36734694
    ,-2.20408163
    ,-2.04081633
    ,-1.87755102
    ,-1.71428571
    ,-1.55102041
    ,-1.3877551
    ,-1.2244898
    ,-1.06122449
    ,-0.89795918
    ,-0.73469388
    ,-0.57142857
    ,-0.40816327
    ,-0.24489796
    ,-8.163265e-2
    ,8.163265e-2
    ,0.24489796
    ,0.40816327
    ,0.57142857
    ,0.73469388
    ,0.89795918
    ,1.06122449
    ,1.2244898
    ,1.3877551
    ,1.55102041
    ,1.71428571
    ,1.87755102
    ,2.04081633
    ,2.20408163
    ,2.36734694
    ,2.53061224
    ,2.69387755
    ,2.85714286
    ,3.02040816
    ,3.18367347
    ,3.34693878
    ,3.51020408
    ,3.67346939
    ,3.83673469
    ,4.0]

xTrain :: Matrix D Double
xTrain = delay $ fromUnboxed (Z :. 50 :. 1) (V.fromList xTrainList)

yTrain :: Matrix D Double
yTrain = smap sin xTrain

testTrainingData :: GPTrainingData Double
testTrainingData = GPTrainingData xTrain yTrain

meanTest :: GPMean Double
meanTest matrix = delay (lkT `mulS` (delay $ (fromMaybeSilly (linearSolveS l yTrain))))
    where
        k = kernelFunction matrix xTrain
        l = cholD . symS $ k
        lK = fromMaybeSilly (linearSolveS l k)
        lkT = transpose lK
        fromMaybeSilly = fromMaybe (error "something went wrong")

testGaussianProcess :: GaussianProcess Double
testGaussianProcess = GaussianProcess kernelFunction meanTest
