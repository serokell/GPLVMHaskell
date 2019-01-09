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
import Data.Array.Repa
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
  :: RandomGen gen
  => Matrix D Double   -- ^ learning vectors
  -> Int               -- ^ desired dimention
  -> Maybe Double      -- ^ the value of stop difference of the likelihood expectation between interations of EM.
  -> Maybe Int         -- ^ number of iterations
  -> gen
  -> PPCA
makePPCA leaningVectors desiredDimention stopValue numberOfIterations generator =
  let _learningData@(ADelayed (Z :. d1 :. n) _) = leaningVectors
      _desiredDimentions = desiredDimention
      _numberOfIterations = _numberOfIterations
      _expectationStopDiff = stopValue
      initMatrix = randomMatrixD generator (desiredDimention,desiredDimention)
      initDispersion = fromIntegral $ fst . next $ generator
      (_W, _variance, _finalExpLikelihood) =
        emSteps _learningData initMatrix initDispersion _expectationStopDiff numberOfIterations
  in PPCA{..}

-- TODO: add comments to the most of following local binds
-- Most of them does not have informative names because they are just intermediate calculations.  

emSteps
  :: Matrix D Double
  -> Matrix D Double
  -> Double
  -> Maybe Double
  -> Maybe Int
  -> (Matrix D Double, Double, Double)
emSteps learnMatrix@(ADelayed (Z :. d1 :. n) _) initMatrix@(ADelayed (Z :. d2 :. d3) _) initVariance stopValue numberOfIterations =
  let stopParam :: Either Int Double
      stopParam = case stopValue of
        Just stopExpectDiff -> Right stopExpectDiff
        Nothing -> case numberOfIterations of
          Just finalIteration -> Left finalIteration
          Nothing -> error "You should provide stop value of the log likelihood expectation\
            \ or number of iterations"
      stepOne :: Array D DIM2 Double -> Double -> Int -> (Matrix D Double, Double, Double)
      stepOne oldW oldVariance iteration =
        let idM = fromFunction (Z :. d2 :. d2) (\(Z :. x1 :. y1) -> if y1 == x1 then 1.0 else 0.0)
            idK = fromFunction (Z :. d2 :. d2) (\(Z :. x1 :. y1) -> if y1 == x1 then oldVariance else 0.0)
            m = delay $ idK +^ ((transpose oldW) `mulS` oldW)
            b = delay $ transpose oldW `mulS`
              ((map (\x -> x/oldVariance) idM) -^
                (map (\x -> x/oldVariance) $ (delay $ oldW `mulS` (delay $ invS m)) `mulS` (transpose oldW)))
            xN i = extend (Any :. All :. (1::Int)) $ slice learnMatrix (Any :. i)
            bXn i = delay $ b `mulS` (xN i)
            sigma = idM -^ (b `mulS` oldW)
            expZtrZ i = sigma +^ ((bXn i) `mulS` (transpose $ bXn i))
        in stepTwo (bXn, expZtrZ, xN) oldW oldVariance iteration
      stepTwo :: (Int -> Matrix D Double, Int -> Matrix D Double, Int -> Matrix D Double) -> Matrix D Double -> Double -> Int -> (Matrix D Double, Double, Double)
      stepTwo (bXn, expZtrZ, xN) oldW oldVariance iteration =
        let sumXnEznZ = transpose $ sumListMatrices [ms | nums <- [1..n], let ms = delay $ (xN nums) `mulS` (bXn nums)]
            sumEznZTrZ = sumListMatrices [ms | nums <- [1..n], let ms = delay $ expZtrZ nums]
            newW = delay $ sumXnEznZ `mulS` (delay $ invS sumEznZTrZ)
            secondDenominatorInVariance i = map (*2) $ ((delay $ (transpose $ bXn i) `mulS` newW)) `mulS` (xN i)
            thirdDenominatorInVariance i = (delay $ (expZtrZ i) `mulS` newW) `mulS` (transpose newW)
            varianceStep i = (flip index) (Z :. 0 :. 0) $ ((xN i) `mulS` (transpose $ xN i)) -^ secondDenominatorInVariance i +^ thirdDenominatorInVariance i
            newVariance = (sum [vss | vss <- U.map varianceStep [1..n]])/(fromIntegral $ n*d2)
            maxDiffNewOldW = foldAllS max 0.0 $ map abs $ oldW -^ newW
            diffVariance = abs $ newVariance - oldVariance
            expLogLikelihood = ((fromIntegral $ n*d2)/2.0)*((-1.0)*(log newVariance) - 1.0)
        in case stopParam of
              Left maxIteration -> if maxIteration > iteration then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood)
              Right stopDiffBetweenIterations -> if (max diffVariance maxDiffNewOldW) > stopDiffBetweenIterations then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood)
   in stepOne initMatrix initVariance 0
