module GPLVM.PPCA
       ( PPCA (..)
       , makePPCA
       ) where

import Prelude (isNaN, log)

import Universum hiding (All, Any, Vector, map, toList, transpose, (++))

import qualified Universum as U

import GPLVM.Types
import GPLVM.Util

import Data.Array.Repa
import Data.Array.Repa.Repr.ForeignPtr
import Data.List ((\\))
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector, diag)
import System.Random

data PPCA = PPCA
  {  _hasMissedData       :: Bool
   , _learningData        :: Matrix D Double
   , desiredDimension   :: Int
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
makePPCA leaningVectors _hasMissedData desiredDimension stopParameter generator =
  let _learningData@(ADelayed (Z :. n :. _) _) = leaningVectors
      initMatrix = randomMatrixD generator (n, desiredDimension)
      initDispersion = fromIntegral $ fst . next $ generator
      (_W, _variance, _finalExpLikelihood, _restoredMatrix) =
        if _hasMissedData
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
  -> (Matrix D Double, Double, Double, Maybe (Matrix D Double))
emStepsFast learnMatrix@(ADelayed (Z :. _ :. n) _) initMatrix@(ADelayed (Z :. d :. _) _) initVariance stopParam =
  let learnMatrixCentered = transpose $ substractMean (transpose learnMatrix)

      stepOne :: (HasCallStack) => Array D DIM2 Double -> Double -> Int -> (Matrix D Double, Double, Double, Maybe (Matrix D Double))
      stepOne (!oldW) oldVariance iteration =
        let m = delay $ mapDiagonal (oldVariance +) $ (transpose oldW) `mulS` oldW
            invM = invS m
            expEznZ = delay $ invM `mul` (transpose oldW `mulS` learnMatrixCentered)
            expXnEZnZ = learnMatrixCentered `mulS` (transpose expEznZ)
            expZtrZ = computeS $ (map (*((fromIntegral n)*oldVariance)) $ delay invM) +^ (expEznZ `mulS` (transpose $ expEznZ)) :: Matrix F Double
        in stepTwo (expEznZ, expZtrZ, expXnEZnZ) oldW oldVariance iteration

      stepTwo :: (HasCallStack) => (Matrix D Double, Matrix F Double, Matrix F Double) -> Matrix D Double -> Double -> Int -> (Matrix D Double, Double, Double, Maybe (Matrix D Double))
      stepTwo (expEznZ, expZtrZ, expXnEZnZ) oldW oldVariance iteration =
        let newW = delay $ expXnEZnZ `mul` (inv expZtrZ)
            u = chol $ trustSym expZtrZ
            wr = newW `mulS` (transpose u)
            totalSum = sumAllS $ map (^ (2 :: Int)) learnMatrixCentered
            secondDenominatorInVariance = (*2) $ sumAllS $ expEznZ *^ (delay $ (transpose newW) `mulS` learnMatrixCentered)
            thirdDenominatorInVariance  = sumAllS $ map (^(2 :: Int)) wr
            newVariance = (totalSum - secondDenominatorInVariance + thirdDenominatorInVariance)/(fromIntegral $ n*d)
            maxDiffNewOldW = foldAllS max 0.0 $ map abs $ oldW -^ newW
            diffVariance = abs $ newVariance - oldVariance
            newC = mapDiagonal ((+) newVariance) $ newW `mulS` (transpose newW)
            newInvM = delay $ invS $ mapDiagonal (newVariance +) $ (transpose newW) `mulS` newW
            newInvC = (mapDiagonal (+(1/newVariance))
                        (map (\x -> -x/newVariance) $ (delay $ newW `mulS` newInvM) `mulS` (transpose newW)))
            expLogLikelihood = (-1.0)*((fromIntegral n)/2.0)*((fromIntegral d)*(log $ 2*pi)
              + (log $ detS newC)
              + (trace2S $ computeS $ map (/(fromIntegral $ n-1)) $ newInvC `mulS` (delay $ learnMatrixCentered `mulS` (transpose learnMatrixCentered))))
        in case stopParam of
              Left maxIteration -> if maxIteration > iteration then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood, Nothing)
              Right stopDiffBetweenIterations -> if (max diffVariance maxDiffNewOldW) > stopDiffBetweenIterations then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood, Nothing)

  in stepOne initMatrix initVariance 0

emStepsMissed
  :: (HasCallStack) => Matrix D Double
  -> Matrix D Double
  -> Double
  -> Either Int Double
  -> (Matrix D Double, Double, Double, Maybe (Matrix D Double))
emStepsMissed learnMatrix@(ADelayed (Z :. _ :. n) _) initMatrix@(ADelayed (Z :. d :. _) _) initVariance stopParam =
  let initMu = delay $ fromListUnboxed (Z:.d:.1) $ replicate d 0.0 -- extend (Any :. (1 :: Int)) $ meanColumnWithNan learnMatrix

      stepOne :: (HasCallStack) => Array D DIM2 Double -> Array D DIM2 Double -> Double -> Int -> (Matrix D Double, Double, Double, Maybe (Matrix D Double))
      stepOne oldMu (!oldW) oldVariance iteration =
        let yi i = extend (Any :. (1 :: Int)) $ slice learnMatrix (Any :. (i :: Int))
            yList i = toList (yi i)
            unknownIndices i =
              let zipped = zip [0..] (yList i)
              in U.map fst $ filter (isNaN . snd) zipped
            yPi i  = deleteRows (unknownIndices i) (yi i)
            muP i =  deleteRows (unknownIndices i) oldMu
            oldWP i = deleteRows (unknownIndices i) oldW
            mP i = mapDiagonal (+ oldVariance) $ (transpose (oldWP i)) `mulS` (oldWP i)
            invMP i = delay $ invS $ mP i
            expXi i = delay $ (delay $ (delay $ invMP i) `mulS` (transpose (oldWP i))) `mulS` (yPi i -^ muP i)
        in stepTwo ([yPi, expXi, invMP], unknownIndices, oldW) oldVariance iteration

      stepTwo :: (HasCallStack) => ([Int -> Matrix D Double], Int -> [Int], Matrix D Double) -> Double -> Int -> (Matrix D Double, Double, Double, Maybe (Matrix D Double))
      stepTwo ([yP, expXi, invMP], unknownIndices, oldW) oldVariance iteration =
        let expX = foldl1 (\acc i -> acc ++ i) $ U.map expXi [0..(n-1)]
            newMu = extend (Any :. (1 :: Int)) $ meanColumnWithNan $ (learnMatrix -^ (oldW `mulS` expX))

            newW =
              let yj j = extend (Any :. (1 :: Int) :. All) $ slice learnMatrix (Any :. (j :: Int) :. All)
                  yJList j = toList (yj j)
                  unknownColumns j =
                    let zipped = zip [0..] (yJList j)
                    in U.map fst $ filter (isNaN . snd) zipped
                  knownColumns j = [0..(n-1)] \\ (unknownColumns j)
                  expXPj j = deleteColumns (unknownColumns j) expX
                  sumInvMP j = sumListMatrices $ U.map invMP (knownColumns j)
                  gj j = ((expXPj j) `mulS` (transpose $ expXPj j)) +^ (map (*oldVariance) (sumInvMP j))
                  yPj j = deleteColumns (unknownColumns j) (yj j)
                  expXtrXj j = (expXPj j) `mulS` (transpose (map (\x -> x - (newMu ! (Z:.j:.0))) (yPj j)))
                  newWj j = (gj j) `solveS` (delay $ expXtrXj j)
              in transpose $ foldl1 (\acc j -> acc ++ j) $ U.map (\x -> delay $ newWj x) [0..(d-1)]

            newVariance =
              let newWPi i = deleteRows (unknownIndices i) newW
                  newMuP i = deleteRows (unknownIndices i) newMu
                  varianceStep i = sumAllS $
                    (map (^(2 :: Int)) (yP i -^ ((newWPi i) `mulS` (expXi i)) -^ (newMuP i)))
                    +^ (map (*oldVariance) $ diag $ delay $ (delay $ (delay $ newWPi i) `mulS` (invMP i)) `mulS` (transpose $ newWPi i))
                  knownVars = foldAllS (\acc x -> if isNaN x then acc else acc + 1.0) 0.0 learnMatrix
              in (sum $ U.map varianceStep [0..(n-1)])/knownVars

            expLogLikelihood =
              let newMuP i = deleteRows (unknownIndices i) newMu
                  newWPi i = deleteRows (unknownIndices i) newW
                  yCenteredP i = (yP i) -^ (newMuP i)
                  yCenteredProd i = delay $ (yCenteredP i) `mulS` (transpose $ yCenteredP i)
                  invMY i = mapDiagonal (+newVariance) $ (newWPi i) `mulS` (transpose $ newWPi i)
                  knownIndices i = [0..(d-1)] \\ (unknownIndices i)
                  expLogLikelihoodStep i = (fromIntegral $ length $ knownIndices i)*(log $ 2*pi) + (log $ detS $ invMY i) +
                    (trace2S $ computeS $ delay $ (invMY i) `solveS` (yCenteredProd i))
              in (-1.0)*(sum $ U.map expLogLikelihoodStep [0..(n-1)])/2.0

            restoredData =
              let muX = extend (Any :. (1 :: Int)) $ meanColumn (transpose expX)
                  finalMu = slice (newMu +^ (newW `mulS` muX)) (Any :. (0 :: Int))
                  wtrw = (transpose newW) `mulS` newW
                  factor1 = newW `mulS` (delay $ inv wtrw)
                  factor2 = computeS $ mapDiagonal (+newVariance) wtrw
                  restoredCentered =  factor1 `mul` factor2 `mul` (computeS $ expX)
              in (delay restoredCentered) +^ (extend (Any :. (n :: Int)) finalMu)

            diffVariance = abs $ newVariance - oldVariance
            maxDiffNewOldW = foldAllS max 0.0 $ map abs $ oldW -^ newW

        in case stopParam of
              Left maxIteration ->
                if maxIteration > iteration
                then stepOne newMu newW newVariance (iteration + 1)
                else (newW,newVariance, expLogLikelihood, Just restoredData)
              Right stopDiffBetweenIterations ->
                if (max diffVariance maxDiffNewOldW) > stopDiffBetweenIterations
                then stepOne newMu newW newVariance (iteration + 1)
                else (newW,newVariance, expLogLikelihood, Just restoredData)
  in stepOne initMu initMatrix initVariance 0

{-
emStepsMissed
  :: (HasCallStack) => Matrix D Double
  -> Matrix D Double
  -> Double
  -> Either Int Double
  -> (Matrix D Double, Double, Double)
emStepsMissed learnMatrix@(ADelayed (Z :. d1 :. n) _) initMatrix@(ADelayed (Z :. d :. d2) _) initVariance stopParam =
  let initMu = extend (Any :. (1 :: Int)) $ meanColumnWithNan learnMatrix

      stepOne :: (HasCallStack) => Array D DIM2 Double -> Double -> Int -> Array D DIM2 Double -> (Matrix D Double, Double, Double)
      stepOne (!oldP) oldVariance iteration oldMu =
        let -- expT :: HasCallStack => Int -> Matrix D Double
            expT i = delay $ (delay $ invW i) `mulS`(transpose $ sumPjXsubMu i)
--            expT i = trace ("expTi: " ++ (show @String $ extent $ expT1 i)) $ expT1 i
            x i = extend (Any :. (1 :: Int)) $ slice learnMatrix (Any :. i)
            xList i = toList (x i)
            knownIndices i =
              let zipped = zip [0..] (xList i)
              in U.map fst $ filter (not . isNaN . snd) zipped
            oldPJ j = extend (Any :. (1 :: Int) :. All) $ slice oldP (Any :. j :. All)
 --           oldPJ j = trace ((show @String j) ++ " hui : " ++ (show @String (computeS $ delay $ oldPJ1 j :: Array U DIM2 Double) ) ) $ oldPJ1 j
            -- oldPP i = getRows (knownIndices i) oldP --oldMu ! (Z :. j :. 0)
            xiP i = getRows (knownIndices i) (x i) -- !  (Z :. j :. 0)
            oldPP :: HasCallStack => Int -> Matrix D Double
            oldPP i = getRows (knownIndices i) oldP
            normOldPP i = sumAllS $ map (^2) $ oldPP i
            oldPPJ :: HasCallStack => Int -> Int -> Matrix D Double
            oldPPJ i j = extend (Any :. (1 :: Int)) $ (slice (oldPP i) (Any :. j) :: Vector D Double)
--            oldPPJ i j = trace ((show @String 89) ++ " : " ++ (show @String (computeS $ delay $ (oldPPJ1 i j) :: Array U DIM2 Double) ) ) $ oldPPJ1 i j
            w :: HasCallStack => Int -> Matrix D Double
            w i = map (+ (normOldPP i)) $ fromFunction (Z:. d :. d) (\(Z :. y :. x) -> if y /= x then 0.0 else oldVariance) --mapDiagonal (+ oldVariance) $ (oldPP i) `mulS` (transpose $ oldPP i)
            invW i = invS $ w i
            invWJ i j = extend (Any :. (1 :: Int)) $ (slice (invW i) (Any :. (j :: Int)) :: Vector D Double )
            sumPjXsubMu :: (HasCallStack) => Int -> Matrix D Double
            sumPjXsubMu i = sumListMatrices $ U.map (\j ->
--              trace (("oldPPJ" :: String) ++ " : " ++ (show @String (computeS $ delay $ (oldPPJ i j) :: Array U DIM2 Double) ) ) $
                map (*((learnMatrix ! (Z:. j :. i)) - (oldMu ! (Z:.j:.0)))) (oldPJ j)) (knownIndices i) -- sumListMatrices $ U.map (\j -> map (* ((xij i j) - (oldMu j))) (oldP j)) knownIndices
--            sumPjXsubMu i = sumPjXsubMu1 i

            expX = fromFunction (Z :. d1 :. n)
                    (\(Z :. j :. i) -> if isNaN $ learnMatrix ! (Z :. j :. i)
                                       then (((oldPJ j) `mulS` (delay $ expT i)) ! (Z :. 0 :. 0)) + (oldMu ! (Z :. j :. 0))
                                       else learnMatrix ! (Z :. j :. i))
            expXi i = extend (Any :. (1 :: Int)) $ slice expX (Any :. i)
            expXij :: Int -> Int -> Double
            expXij i j = expX ! (Z :. j :. i) -- if isNaN $ learnMatrix ! (Z :. j :. i) then (((oldPJ j) `mulS` (delay $ expT i)) ! (Z :. 0 :. 0)) + (oldMu ! (Z :. j :. 0)) else learnMatrix ! (Z :. j :. i)
--            expXi i = fromFunction (Z :. d1 :. 1) (\(Z :. y :. x) -> expXij i y)

            tiTrTi i = map (* oldVariance) (invW i) +^ ((delay $ expT i) `mulS` (transpose $ expT i))

            expXiTrXi i = fromFunction (Z :. d :. d)
              (\(Z :. j :. k) ->
                if | (isNaN $ learnMatrix ! (Z :. j :. i)) && (isNaN $ learnMatrix ! (Z :. k :. i)) && j /= k ->
                       (((delay $ (oldPJ j) `mulS` (delay $ invW i)) `mulS` (transpose $ oldPJ k)) ! (Z :. 0 :. 0))*oldVariance + (expXij i j)*(expXij i k)
                   | (isNaN $ learnMatrix ! (Z :. j :. i)) && (isNaN $ learnMatrix ! (Z :. k :. i)) && j == k ->
                        (((delay $ (oldPJ j) `mulS` (delay $ invW i)) `mulS` (transpose $ oldPJ k)) ! (Z :. 0 :. 0) + 1.0)*oldVariance + (expXij i j)*(expXij i k)
                   | (isNaN $ learnMatrix ! (Z :. j :. i)) && (not $ isNaN $ learnMatrix ! (Z :. k :. i)) ->
                        (expXij i j)*(learnMatrix ! (Z :. k :. i))
                   | (not $ isNaN $ learnMatrix ! (Z :. j :. i)) && (isNaN $ learnMatrix ! (Z :. k :. i)) ->
                        (learnMatrix ! (Z :. j :. i))*(expXij i k)
                   | (not $ isNaN $ learnMatrix ! (Z :. j :. i)) && (not $ isNaN $ learnMatrix ! (Z :. k :. i)) ->
                        (learnMatrix ! (Z :. j :. i))*(learnMatrix ! (Z :. k :. i)) )
            expXiTrTi i = fromFunction (Z :. d :. d)
              (\(Z :. k :. j) ->
                if   (isNaN $ learnMatrix ! (Z :. j :. i))
                then (((oldP) `mulS` (delay $ invW i)) ! (Z :. 0 :. j))*(oldVariance) +
                  (((expXi i) `mulS` (transpose $ expT i)) ! (Z :. 0 :. 0))
                else (x i ! (Z:. j :. 0))*((expT i) ! (Z :. k :. 0)))

        in stepTwo ([x, expT, oldPJ, expXi, tiTrTi, expXiTrXi, expXiTrTi], knownIndices, oldP) iteration oldVariance

      stepTwo :: (HasCallStack) => ([Int -> Matrix D Double], Int -> [Int], Matrix D Double) -> Int -> Double -> (Matrix D Double, Double, Double)
      stepTwo ([x, expT, oldPJ, expXi, tiTrTi, expXiTrXi, expXiTrTi], knownIndices, oldP) iteration oldVariance =
        let newMu :: Matrix D Double
            newMu = map (\x -> x/(fromIntegral n)) $ sumListMatrices $ U.map (\i -> (expXi i) -^ (oldP `mulS` (expT i))) [0..(n - 2)]
--            newMu :: Matrix D Double
--            newMu = trace ("nuewMu" ++ (show @String (computeS newMu1 :: Array U DIM2 Double))) newMu1
            newP = delay $ (sumListMatrices $ U.map (\i -> (expXiTrTi i) -^ (newMu `mulS` (transpose $ expT i))) [0..(n-2)])
                     `mulS` (delay $ invS $ sumListMatrices $ U.map (\i -> tiTrTi i) [0..(n-2)])
            varianceStep i = map (+ ((newMu `mulS` (transpose newMu) ! (Z :. 0 :. 0)) - 2*((newMu `mulS` (transpose $ expXi i)) ! (Z :. 0 :. 0)) ) ) $
              (expXiTrXi i) -^ (map (*2) $ (expXiTrTi i) `mulS` (transpose $ newP))
                +^ (map (*(2 * (newMu `mulS` (transpose $ expT i)) ! (Z:.0:.0))) (transpose newP)) +^ (delay (newP `mulS` (tiTrTi i)) `mulS` (transpose $ newP))
            newVariance = (sum $ U.map (trace2S . computeS . varianceStep) [0..(n-1)])/(fromIntegral $ n*d)
            diffVariance = abs $ newVariance - oldVariance
            maxDiffNewOldW = foldAllS max 0.0 $ map abs $ oldP -^ newP
            expLogLikelihood = 1.0
        in case stopParam of
              Left maxIteration -> if maxIteration > iteration then stepOne newP newVariance (iteration + 1) newMu else (newP,newVariance, expLogLikelihood)
              Right stopDiffBetweenIterations -> if (max diffVariance maxDiffNewOldW) > stopDiffBetweenIterations then stepOne newP newVariance (iteration + 1) newMu else (newP,newVariance, expLogLikelihood)
  in stepOne initMatrix initVariance 0 initMu -}
