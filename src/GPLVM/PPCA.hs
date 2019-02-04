module GPLVM.PPCA
       ( PPCA (..)
       , makePPCA
       ) where

import Prelude (isNaN, log)

import Universum hiding (All, Any, Vector, map, toList, transpose)

import qualified Universum as U

import GPLVM.Types
import GPLVM.Util

import Data.Array.Repa hiding ((++))
import Data.Array.Repa.Algorithms.Matrix
import Data.Array.Repa.Repr.ForeignPtr
import Data.List (elemIndex)
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import System.Random

data PPCA = PPCA
  {  _hasMissedData       :: Bool
   , _learningData        :: Matrix D Double
   , desiredDimentions   :: Int
   , stopParameter       :: Either Int Double
   , _variance            :: Double
   , _W                   :: Matrix D Double
   , _finalExpLikelihood  :: Double
   , _restoredMatrix      :: Maybe (Matrix D Double)
   }

makePPCA
  :: (HasCallStack) => RandomGen gen
  => Matrix D Double   -- ^ learning vectors
  -> Bool              -- ^ Is there missed data
  -> Int               -- ^ desired dimention
  -> Either Int Double -- ^ number of iterations or
                       --   the value of stop difference of the likelihood expectation between interations of EM.
  -> gen
  -> PPCA
makePPCA leaningVectors missedData desiredDimentions stopParameter generator =
  let _learningData@(ADelayed (Z :. n :. d1) _) = leaningVectors
--      initMatrix = delay $ fromListUnboxed (Z :. n :. desiredDimention) [-0.925848419891795, 0.181032149310173, 0.0187743520609800, -0.319885618861164, 0.00532168979485395, 0.657479738708321, -0.425960758296183, 0.825430464554896, 1.13914949095078, 0.584245605925614, -2.03402726516440, -0.229292809268995, 0.427361247541929, -1.61645129547965, -1.32855255183639, 0.262250813294937] -- randomMatrixD generator (n, desiredDimention)
      initMatrix = randomMatrixD generator (n, desiredDimentions)
      initDispersion = fromIntegral $ fst . next $ generator -- 883.762
      (_W, _variance, _finalExpLikelihood) =
        if missedData
        then emStepsMissed _learningData initMatrix initDispersion stopParameter
        else emStepsFast _learningData initMatrix initDispersion stopParameter
  in PPCA{..}

-- TODO: add comments to the most of following local binds
-- Most of them does not have informative names because they are just intermediate calculations.

emStepsFast
  :: (HasCallStack) => Matrix D Double
  -> Matrix D Double
  -> Double
  -> Either Int Double
  -> (Matrix D Double, Double, Double)
emStepsFast learnMatrix@(ADelayed (Z :. d1 :. n) _) initMatrix@(ADelayed (Z :. d3 :. d2) _) initVariance stopParam =
  let learnMatrixCentered = transpose $ substractMean (transpose learnMatrix)

      stepOne :: (HasCallStack) => Array D DIM2 Double -> Double -> Int -> (Matrix D Double, Double, Double)
      stepOne (!oldW) oldVariance iteration =
        let m = delay $ mapDiagonal (oldVariance +) $ (transpose oldW) `mulS` oldW
            invM = invS m
            expEznZ = delay $ invM `mul` (transpose oldW `mulS` learnMatrixCentered)
            expXnEZnZ = learnMatrixCentered `mulS` (transpose expEznZ)
            expZtrZ = computeS $ (map (*((fromIntegral n)*oldVariance)) $ delay invM) +^ (expEznZ `mulS` (transpose $ expEznZ)) :: Matrix F Double
        in stepTwo (expEznZ, expZtrZ, expXnEZnZ) oldW oldVariance iteration

      stepTwo :: (HasCallStack) => (Matrix D Double, Matrix F Double, Matrix F Double) -> Matrix D Double -> Double -> Int -> (Matrix D Double, Double, Double)
      stepTwo (expEznZ, expZtrZ, expXnEZnZ) oldW oldVariance iteration =
        let newW = delay $ expXnEZnZ `mul` (inv expZtrZ)
            u = chol $ trustSym expZtrZ
            wr = newW `mulS` (transpose u)
            totalSum = sumAllS $ map (^2) learnMatrixCentered
            secondDenominatorInVariance = (*2) $ sumAllS $ expEznZ *^ (delay $ (transpose newW) `mulS` learnMatrixCentered)
            thirdDenominatorInVariance  = sumAllS $ map (^2) wr
            newVariance = (totalSum - secondDenominatorInVariance + thirdDenominatorInVariance)/(fromIntegral $ n*d3)
            maxDiffNewOldW = foldAllS max 0.0 $ map abs $ oldW -^ newW
            diffVariance = abs $ newVariance - oldVariance
            newC = mapDiagonal ((+) newVariance) $ newW `mulS` (transpose newW)
            newInvM = delay $ invS $ mapDiagonal (newVariance +) $ (transpose newW) `mulS` newW
            newInvC = (mapDiagonal (+(1/newVariance))
                        (map (\x -> -x/newVariance) $ (delay $ newW `mulS` newInvM) `mulS` (transpose newW)))
            expLogLikelihood = (-1.0)*((fromIntegral n)/2.0)*((fromIntegral d3)*(log $ 2*pi)
              + (log $ detS newC)
              + (trace2S $ computeS $ map (/(fromIntegral $ n-1)) $ newInvC `mulS` (delay $ learnMatrixCentered `mulS` (transpose learnMatrixCentered))))
        in case stopParam of
              Left maxIteration -> if maxIteration > iteration then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood)
              Right stopDiffBetweenIterations -> if (max diffVariance maxDiffNewOldW) > stopDiffBetweenIterations then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood)

   in stepOne initMatrix initVariance 0

emStepsMissed
  :: (HasCallStack) => Matrix D Double
  -> Matrix D Double
  -> Double
  -> Either Int Double
  -> (Matrix D Double, Double, Double)
emStepsMissed learnMatrix@(ADelayed (Z :. d1 :. n) _) initMatrix@(ADelayed (Z :. d3 :. d2) _) initVariance stopParam =
  let learnMatrixCentered = transpose $ substractMeanWithNan (transpose learnMatrix) :: Array D DIM2 Double

      stepOne :: (HasCallStack) => Array D DIM2 Double -> Double -> Int -> (Matrix D Double, Double, Double)
      stepOne (!oldW) oldVariance iteration =
        let y i = extend (Any :. (1 :: Int)) $ slice learnMatrixCentered (Any :. (i :: Int))
            yList i = toList (y i)
            unknownIndices i =
              let zipped = zip [0..] (yList i)
              in U.map fst $ filter (isNaN . snd) zipped
            yP i = deleteRows (unknownIndices i) (y i)
            yM i = getRows (unknownIndices i) (y i)
            oldWP i = deleteRows (unknownIndices i) oldW
            oldWM i = getRows (unknownIndices i) oldW
            mP i = mapDiagonal (+ oldVariance) $ (transpose (oldWP i)) `mulS` (oldWP i)
            invMP i = delay $ invS $ mP i
            expX i = delay $ (delay $ (delay $ invMP i) `mulS` (transpose (oldWP i))) `mulS` (yP i)
            expYM i = delay $ (oldWM i) `mulS` (expX i)
            expXiTrXi i = (map (*oldVariance) (invMP i)) +^ ((expX i) `mulS` (transpose $ expX i))
            expYMiTrYMi i = (map (*oldVariance) $ (mapDiagonal (+1) $ delay $ (delay $ (oldWM i) `mulS` (invMP i)) `mulS` (transpose $ oldWM i)))
              +^ ((expYM i) `mulS` (transpose $ expYM i))
            expYMiXi i = (map (*((-1)*oldVariance)) $ (oldWM i) `mulS` (invMP i)) +^ (expYM i) `mulS` (transpose $ expX i)
        in stepTwo ([yP, yM, oldWP, expX, expYM, expXiTrXi, expYMiTrYMi, expYMiXi, y], unknownIndices, oldW) oldVariance iteration

      stepTwo :: (HasCallStack) => ([Int -> Matrix D Double], Int -> [Int], Matrix D Double) -> Double -> Int -> (Matrix D Double, Double, Double)
      stepTwo ([yP, yM, oldWP, expX, expYM, expXiTrXi, expYMiTrYMi, expYMiXi, learnColumn], unknownIndices, oldW) oldVariance iteration =
        let g i = fromFunction (Z :. d1 :. d3) (\(Z :. x :. y) ->
                    if (elem x (unknownIndices i))
                    then ((expYMiXi i) ! (Z :. (getUnknownIndex i x) :. y))
                    else ((learnColumn i) ! (Z :. x :. 0))*((expX i) ! (Z :. y :. 0)) )
            trcs i = "expYM i: " ++ (show @String (computeS $ (expYM i) :: Array U DIM2 Double)) ++ "\n" ++ "newWP i: "
                        ++ (show @String (computeS $ (newWP i) :: Array U DIM2 Double)) ++ "\n" ++ "expX i: "
                        ++ (show @String (computeS $ (expX i) :: Array U DIM2 Double)) ++ "\n" ++ "expYMiXi i: "
                        ++ (show @String (computeS $ (expYMiXi i) :: Array U DIM2 Double)) ++ "\n" ++ "newWM i: "
                        ++ (show @String (computeS $ (newWM i) :: Array U DIM2 Double)) ++ "\n" ++ "lalala: "
                        ++ (show @String $ (if (length $ toList $ expYMiXi i) /= 0
                                              then ((*2) $ trace2S $ computeS $ delay $ (transpose $ newWM i) `mulS` (expYMiXi i))
                                              else 0.0))
            getUnknownIndex i x = fromMaybe (error "No such element in unknown observations") $ elemIndex x (unknownIndices i)
            columnIndices = [0..(n-1)]
            sumOfGi = sumListMatrices $ U.map g columnIndices
            sumOfCovXInv = delay $ invS $ sumListMatrices $ U.map expXiTrXi columnIndices
            newW = delay $ sumOfGi `mulS` sumOfCovXInv
            newWP i = deleteRows (unknownIndices i) newW
            newWM i = getRows (unknownIndices i) newW
            varianceStep i = ((transpose $ yP i) `mulS` (yP i) ! (Z :. 0 :. 0))
                               - ((*2) $ ((delay $ (transpose $ expX i) `mulS` (transpose $ newWP i)) `mulS` (yP i)) ! (Z :. 0 :. 0))
                               - (if (length $ toList $ expYMiXi i) /= 0
                                  then ((*2) $ trace2S $ computeS $ delay $ (transpose $ newWM i) `mulS` (expYMiXi i))
                                  else 0.0)
                               + (trace2S $ computeS $ delay $ (delay ((transpose newW) `mulS` newW)) `mulS` (expXiTrXi i))
                               + (trace2S $ computeS $ expYMiTrYMi i)
            newVariance1 = (sum $ U.map varianceStep columnIndices)/(fromIntegral $ n * d3)
            newVariance = trace ("newVariance: " ++ (show @String newVariance1) ++ "\n" ++ "yUnknown: " ++ (show @String (computeS $ expYM 2 :: Array U DIM2 Double))) newVariance1
            diffVariance = trace ("23" :: String) $ abs $ newVariance - oldVariance
            maxDiffNewOldW = trace ("24" :: String) $ foldAllS max 0.0 $ map abs $ oldW -^ newW
            expLogLikelihood = 1.0
        in case stopParam of
              Left maxIteration -> if maxIteration > iteration then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood)
              Right stopDiffBetweenIterations -> if (max diffVariance maxDiffNewOldW) > stopDiffBetweenIterations then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood)
  in stepOne initMatrix initVariance 0
