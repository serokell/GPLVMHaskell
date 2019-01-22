module GPLVM.PPCA
       ( PPCA (..)
       , makePPCA
       ) where

import Prelude (log)
import Universum hiding (All, Any, Vector, map, transpose)
import qualified Universum as U (map)

import GPLVM.Types
import GPLVM.Util

import Control.Lens (makeLenses)
import Data.Array.Repa hiding ((++))
import Data.Array.Repa.Algorithms.Matrix
import GHC.Stack
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import System.Random

data PPCA = PPCA
  {  _learningData        :: Matrix D Double
   , _desiredDimentions   :: Int
   , _numberOfIterations  :: Maybe Int
   , _expectationStopDiff :: Maybe Double
   , _variance            :: Double
   , _W                   :: Matrix D Double
   , _finalExpLikelihood  :: Double
   }

makePPCA
  :: (HasCallStack) => RandomGen gen
  => Matrix D Double   -- ^ learning vectors
  -> Int               -- ^ desired dimention
  -> Maybe Double      -- ^ the value of stop difference of the likelihood expectation between interations of EM.
  -> Maybe Int         -- ^ number of iterations
  -> gen
  -> PPCA
makePPCA leaningVectors desiredDimention stopValue numberOfIterations generator =
  let _learningData@(ADelayed (Z :. n :. d1) _) = leaningVectors
      _desiredDimentions = desiredDimention
      _numberOfIterations = _numberOfIterations
      _expectationStopDiff = stopValue
      initMatrix = delay $ fromListUnboxed (Z :. n :. desiredDimention) [-1.155170597359062, -0.16620412451513678, 0.25738948150146834, -0.14530511413838565, 1.2142412622490049, 0.8594441985274117, 0.2533011836594896, 1.0258331618344974, 0.2985607079847031, -0.6557949684137585, 0.2863661276273283, -0.8279110126834314, 9.479438827868428e-2, -0.24701919460237498, 1.6684212084781769, -0.8866102872137365] --randomMatrixD generator (n, desiredDimention)
      initDispersion = 2.042050972e9 --fromIntegral $ fst . next $ generator
      (_W, _variance, _finalExpLikelihood) =
        emSteps _learningData initMatrix initDispersion _expectationStopDiff numberOfIterations
  in PPCA{..}

-- TODO: add comments to the most of following local binds
-- Most of them does not have informative names because they are just intermediate calculations.

emSteps
  :: (HasCallStack) => Matrix D Double
  -> Matrix D Double
  -> Double
  -> Maybe Double
  -> Maybe Int
  -> (Matrix D Double, Double, Double)
emSteps learnMatrix@(ADelayed (Z :. d1 :. n) _) initMatrix@(ADelayed (Z :. d3 :. d2) _) initVariance stopValue numberOfIterations =
  let stopParam :: Either Int Double
      stopParam = case stopValue of
        Just stopExpectDiff -> Right stopExpectDiff
        Nothing -> case numberOfIterations of
          Just finalIteration -> Left finalIteration
          Nothing -> error "You should provide stop value of the log likelihood expectation\
            \ or number of iterations"
      stepOne :: (HasCallStack) => Array D DIM2 Double -> Double -> Int -> (Matrix D Double, Double, Double)
      stepOne oldW oldVariance iteration =
        let learnMatrixCentered = trace (show @String (computeS learnMatrix :: Array U DIM2 Double)) $ transpose $ substractMean (transpose learnMatrix)
            m1 = delay $ mapDiagonal (oldVariance +) $ (transpose oldW) `mulS` oldW
            m = if iteration == 0 then trace (show @String (computeS m1 :: Array U DIM2 Double)) m1 else m1
            b1 = delay $ transpose oldW `mulS`
                  (mapDiagonal ((1/oldVariance) +)
                    (map (\x -> -x/oldVariance) $ (delay $ oldW `mulS` (delay $ invS m)) `mulS` (transpose oldW)))
            b = if iteration == 0 then trace ("b: " ++ (show @String (computeS b1 :: Array U DIM2 Double))) b1 else b1
            xN1 i = extend (Any :. All :. (1::Int)) $ slice learnMatrixCentered (Any :. i)
            xN i = if i == 0 && iteration == 0
                   then trace ("xN " ++ (show i) ++ " " ++ (show @String $ (computeS $ xN1 i :: Array U DIM2 Double))) $ xN1 i
                   else xN1 i
            bXn1 i = delay $ b `mulS` (xN i)
            bXn i = if i == 0 && iteration == 0
                    then trace ("bxN: " ++ (show @String $ (computeS $ bXn1 i :: Array U DIM2 Double))) (bXn1 i)
                    else (bXn1 i)
            sigma1 = mapDiagonal (\x -> x + 1) $ map (\x -> -x) $ b `mulS` oldW
            sigma = if iteration == 0 then trace ("sigma: " ++ (show @String $ (computeS $ sigma1 :: Array U DIM2 Double))) sigma1 else sigma1
            expZtrZ i = sigma +^ ((bXn i) `mulS` (transpose $ bXn i))
        in stepTwo (bXn, expZtrZ, xN) oldW oldVariance iteration
      stepTwo :: (HasCallStack) => (Int -> Matrix D Double, Int -> Matrix D Double, Int -> Matrix D Double) -> Matrix D Double -> Double -> Int -> (Matrix D Double, Double, Double)
      stepTwo (bXn, expZtrZ, xN) oldW oldVariance iteration =
        let sumXnEznZ = sumListMatrices [ms | nums <- [0..(d1-1)], let ms = delay $ (xN nums) `mulS` (transpose $ bXn nums)]
--            sumXnEznZ = trace ("sumXnEznZ: " ++ (showShape $ extent sumXnEznZ1)) sumXnEznZ1
            sumEznZTrZ = sumListMatrices [ms | nums <- [0..(d1-1)], let ms = delay $ expZtrZ nums]
--            sumEznZTrZ = trace ("sumEznZTrZ " ++ (show @String $ extent sumEznZTrZ1)) sumEznZTrZ1
            newW = delay $ sumXnEznZ `mulS` (delay $ invS sumEznZTrZ)
--            newW = trace ("newW: " ++ (showShape $ extent newW1)) newW1
            secondDenominatorInVariance i = map (*2) $ ((delay $ (transpose $ bXn i) `mulS` (transpose newW))) `mulS` (xN i)
 --           secondDenominatorInVariance i = trace ("sec " ++ (show @String $ (extent $ secondDenominatorInVariance1 i ))) secondDenominatorInVariance1 i
            thirdDenominatorInVariance i = trace2S . computeS . delay $ (delay $ (expZtrZ i) `mulS` (transpose newW)) `mulS` newW
       --     thirdDenominatorInVariance i = trace ("third " ++ (show @String i) ++ " " ++ (show @String $ thirdDenominatorInVariance1 i)) $ thirdDenominatorInVariance1 i
            varianceStep i = thirdDenominatorInVariance i + ((flip index) (Z :. 0 :. 0) $ ((transpose $ xN i) `mulS` (xN i)) -^ secondDenominatorInVariance i)
            newVariance = trace ((show @String (d3*d1)) ++ " : " ++ (show @String $ (sum [vss | vss <- U.map varianceStep [0..(d1-1)]]))) (sum [vss | vss <- U.map varianceStep [0..(d1-1)]])/(fromIntegral $ d1*d3)
            maxDiffNewOldW = foldAllS max 0.0 $ map abs $ oldW -^ newW
            diffVariance = abs $ newVariance - oldVariance
            expLogLikelihood = ((fromIntegral $ d1*d3)/2.0)*((-1.0)*(log newVariance) - 1.0)
        in case stopParam of
              Left maxIteration -> if maxIteration > iteration then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood)
              Right stopDiffBetweenIterations -> if (max diffVariance maxDiffNewOldW) > stopDiffBetweenIterations then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood)
   in stepOne initMatrix initVariance 0
