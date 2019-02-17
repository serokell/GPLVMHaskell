module Math.TestDimensions where

import Data.Array.Repa hiding ((++))
import GHC.TypeLits hiding (someNatVal)
import GPLVM.Types
import GPLVM.Util hiding (trace2S)
import Numeric.Dimensions

import Prelude (log)

import Universum hiding (All, Any, Vector, map, natVal, toList, transpose)

import Data.Array.Repa.Algorithms.Matrix
import Data.Vinyl.TypeLevel (AllConstrained)

import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import System.Random
import Unsafe.Coerce


--TypeFamily MulMatr Dim

newtype DimMatrix r (y :: Nat) (x :: Nat) a
  = DimMatrix { getInternal :: Matrix r a}

data PPCA = PPCA
  {  _learningData        :: Matrix D Double
   , desiredDimentions   :: Int
   , stopParameter       :: Either Int Double
   , _variance            :: Double
   , _W                   :: Matrix D Double
   , _finalExpLikelihood  :: Double
   }

withMat :: Matrix D Double -> (forall x y. (KnownNat x, KnownNat y) => DimMatrix D x y Double -> k) -> k
withMat m f =
    let (Z :. x :. y) = extent m
    in
    case someNatVal (fromIntegral x) of
      SomeNat (Proxy :: Proxy m) -> case someNatVal (fromIntegral y) of
        SomeNat (Proxy :: Proxy n) -> f (DimMatrix @_ @m @n m)

makePPCATypeSafe
  :: RandomGen gen
  => Matrix D Double   -- ^ learning vectors
  -> Int               -- ^ desired dimension
  -> Either Int Double -- ^ number of iterations or
                       --   the value of stop difference of the likelihood expectation between interations of EM.
  -> (forall d x1 y1 x2 y2 . (AllConstrained KnownNat [d, x1, y1, x2, y2], x2 ~ d, y1 ~ y2)
      => DimMatrix D y1 x1 Double
      -> DimMatrix D y2 x2 Double
      -> Double
      -> Either Int Double
      -> (DimMatrix D y2 x2 Double, Double, Double))
  -> gen
  -> PPCA
makePPCATypeSafe leaningVectors desiredDimentions stopParameter func generator =
  let _learningData@(ADelayed (Z :. n :. _) _) = leaningVectors
      initMatrix = randomMatrixD generator (n, desiredDimentions) :: Matrix D Double
      initVariance = fromIntegral $ fst . next $ generator :: Double
      (_W, _variance, _finalExpLikelihood) =
        withMat _learningData $ \(ld :: DimMatrix D y1 x1 Double) ->
        withMat initMatrix $ \(initM :: DimMatrix D y2 x2 Double) ->
        case someNatVal $ fromIntegral desiredDimentions of
          SomeNat (Proxy :: Proxy d) ->
             withEvidence (inferPPCAInputMatrices @d @y1 @y2 @x2) $ convertPPCATypeSafeData $ func @d ld initM initVariance stopParameter
  in PPCA{..}

inferPPCAInputMatrices
  :: forall d y1 y2 x2.
  AllConstrained KnownNat [d, y1, y2, x2]
  => Evidence (x2 ~ d, y1 ~ y2)
inferPPCAInputMatrices
  | natVal (Proxy :: Proxy y1) /= natVal (Proxy :: Proxy y2) =
    error $ toText $ "dimentions y1 and y2 should be equal, but y1 = " ++ (show (natVal (Proxy :: Proxy y1))) ++ " and y2 = " ++ (show (natVal (Proxy :: Proxy y2)))
  | natVal (Proxy :: Proxy x2) /= natVal (Proxy :: Proxy d) =
    error $ toText $ "dimentions x2 and d should be equal, but x2 = " ++ (show (natVal (Proxy :: Proxy x2))) ++ " and d = " ++ (show (natVal (Proxy :: Proxy d)))
  | otherwise = unsafeCoerce (E @((y1 ~ y1), (y1 ~ y1)))

mulM
  :: forall y1 x1 y2 x2.
  ( AllConstrained KnownNat [x1, x2, y1, y2]
  , x1 ~ y2
  )
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> DimMatrix D y1 x2 Double
mulM (DimMatrix m1) (DimMatrix m2) = DimMatrix $ delay $ m1 `mulS` m2

transposeM
  :: (KnownNat y, KnownNat x)
  => DimMatrix D y x Double
  -> DimMatrix D x y Double
transposeM (DimMatrix m) = DimMatrix $ transpose m

mapMM :: (KnownNat y, KnownNat x)
  => (Double -> Double)
  -> DimMatrix D y x Double
  -> DimMatrix D y x Double
mapMM f (DimMatrix m) =  DimMatrix $ map f m

mapDiagonalM :: (KnownNat y, KnownNat x)
  => (Double -> Double)
  -> DimMatrix D y x Double
  -> DimMatrix D y x Double
mapDiagonalM f (DimMatrix m) = DimMatrix $ mapDiagonal f m

invSM
  :: (KnownNat y, KnownNat x, y ~ x)
  => DimMatrix D y x Double
  -> DimMatrix D y x Double
invSM (DimMatrix m) = DimMatrix $ delay $ invS m

substractMeanM
  :: (KnownNat y, KnownNat x)
  => DimMatrix D y x Double
  -> DimMatrix D y x Double
substractMeanM (DimMatrix m) = DimMatrix $ substractMean m

(+^^) :: forall y1 x1 y2 x2.
  ( AllConstrained KnownNat [x1, x2, y1, y2]
  , x1 ~ x2
  , y1 ~ y2
  )
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> DimMatrix D y2 x2 Double
(+^^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 +^ m2

(-^^) :: forall y1 x1 y2 x2.
  ( AllConstrained KnownNat [x1, x2, y1, y2]
  , x1 ~ x2
  , y1 ~ y2
  )
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> DimMatrix D y2 x2 Double
(-^^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 -^ m2

(*^^) :: forall y1 x1 y2 x2.
  ( AllConstrained KnownNat [x1, x2, y1, y2]
  , x1 ~ x2
  , y1 ~ y2
  )
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> DimMatrix D y2 x2 Double
(*^^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 *^ m2

cholM
  :: (KnownNat y, KnownNat x, y ~ x)
  => DimMatrix D y x Double
  -> DimMatrix D y x Double
cholM (DimMatrix m) = DimMatrix $ delay $ chol $ trustSym $ computeS m

sumAllSM
  :: (KnownNat y, KnownNat x)
  => DimMatrix D y x Double
  -> Double
sumAllSM (DimMatrix m) = sumAllS m

foldAllSM :: (KnownNat y, KnownNat x)
  => (Double -> Double -> Double)
  -> Double
  -> DimMatrix D y x Double
  -> Double
foldAllSM f initValue (DimMatrix m) = foldAllS f initValue m

detSM :: (KnownNat y, KnownNat x)
  => DimMatrix D y x Double
  -> Double
detSM (DimMatrix m) = detS m

trace2SM :: (KnownNat y, KnownNat x)
  => DimMatrix D y x Double
  -> Double
trace2SM (DimMatrix m) = trace2S $ computeS m

emStepsFast
  :: forall d y1 x1 y2 x2.
  ( HasCallStack
  , AllConstrained KnownNat [x1, x2, y1, y2]
  , x2 ~ d
  , y1 ~ y2
  )
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> Double
  -> Either Int Double
  -> (DimMatrix D y2 x2 Double, Double, Double)
emStepsFast learnMatrix initMatrix initVariance stopParam =
  let learnMatrixCentered = transposeM $ substractMeanM (transposeM learnMatrix)
      n = natVal (Proxy :: Proxy x1)
      d3 = natVal (Proxy :: Proxy y2)

      stepOne :: (HasCallStack)
              => DimMatrix D y2 x2 Double
              -> Double
              -> Int
              -> (DimMatrix D y2 x2 Double, Double, Double)
      stepOne (!oldW) oldVariance iteration =
        let m = mapDiagonalM (oldVariance +) $ (transposeM oldW) `mulM` oldW
            invM = invSM m
            expEznZ = invM `mulM` (transposeM oldW `mulM` learnMatrixCentered)
            expXnEZnZ = learnMatrixCentered `mulM` (transposeM expEznZ)
            expZtrZ = (mapMM (*((fromIntegral n)*oldVariance)) invM) +^^ (expEznZ `mulM` (transposeM $ expEznZ))
        in stepTwo (expEznZ, expZtrZ, expXnEZnZ) oldW oldVariance iteration

      stepTwo
        :: forall y12 x12 y22 x22 y32 x32 y42 x42 y52 x52 .
           (HasCallStack, x42 ~ d, y12 ~ d, x22 ~ d, y22 ~ d, x32 ~ d, x12 ~ x1, y32 ~ y2, x52 ~ d, y2 ~ y52, y42 ~ y52)
           => (DimMatrix D y12 x12 Double, DimMatrix D y22 x22 Double, DimMatrix D y32 x32 Double)
           -> DimMatrix D y42 x42 Double
           -> Double
           -> Int
           -> (DimMatrix D y52 x52 Double, Double, Double)
      stepTwo (expEznZ, expZtrZ, expXnEZnZ) oldW oldVariance iteration =
        let newW = expXnEZnZ `mulM` (invSM expZtrZ)
            u = cholM expZtrZ
            wr = newW `mulM` (transposeM u)
            totalSum = sumAllSM $ mapMM (^(2 :: Int)) learnMatrixCentered
            secondDenominatorInVariance = (*2.0) $ sumAllSM $ expEznZ *^^ ((transposeM newW) `mulM` learnMatrixCentered)
            thirdDenominatorInVariance  = sumAllSM $ mapMM (^(2 :: Int)) wr
            newVariance = (totalSum - secondDenominatorInVariance + thirdDenominatorInVariance)/(fromIntegral $ n*d3)
            maxDiffNewOldW = foldAllSM max 0.0 $ mapMM abs $ oldW -^^ newW
            diffVariance = abs $ newVariance - oldVariance
            newC = mapDiagonalM ((+) newVariance) $ newW `mulM` (transposeM newW)
            newInvM = invSM $ mapDiagonalM (newVariance +) $ (transposeM newW) `mulM` newW
            newInvC = (mapDiagonalM (+(1/newVariance))
                        (mapMM (\x -> -x/newVariance) $ ((newW `mulM` newInvM) `mulM` (transposeM newW))))
            expLogLikelihood = (-1.0)*((fromIntegral n)/2.0)*((fromIntegral d3)*(log $ 2*pi)
              + (log $ detSM newC)
              + (trace2SM $ mapMM (/(fromIntegral $ n-1)) $ newInvC `mulM` (learnMatrixCentered `mulM` (transposeM learnMatrixCentered))))
        in case stopParam of
              Left maxIteration -> if maxIteration > iteration then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood)
              Right stopDiffBetweenIterations -> if (max diffVariance maxDiffNewOldW) > stopDiffBetweenIterations then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood)

   in stepOne initMatrix initVariance 0
{-
emStepsMissed
  :: (HasCallStack) => Matrix D Double
  -> Matrix D Double
  -> Double
  -> Either Int Double
  -> (Matrix D Double, Double, Double)
emStepsMissed learnMatrix@(ADelayed (Z :. d1 :. n) _) initMatrix@(ADelayed (Z :. d3 :. d2) _) initVariance stopParam =
  let initMu = extend (Any :. (1 :: Int)) $ meanColumnWithNan learnMatrix

      stepOne :: (HasCallStack) => Array D DIM2 Double -> Double -> Int -> Array D DIM2 Double -> (Matrix D Double, Double, Double)
      stepOne (!oldP) oldVariance iteration oldMu =
        let expT :: HasCallStack => Int -> Matrix D Double
            expT i = trace ((show @String i) ++ "kuku : " ++ (show @String (computeS $ delay $ invW i :: Array U DIM2 Double) ) ) $
              delay $ (delay $ invW i) `mulS` (sumPjXsubMu i)
            x i = extend (Any :. (1 :: Int)) $ slice learnMatrix (Any :. i)
            xList i = toList (x i)
            knownIndices i =
              let zipped = zip [0..] (xList i)
              in U.map fst $ filter (not . isNaN . snd) zipped
            oldPJ1 j =  extend (Any :. (1 :: Int)) $ slice oldP (Any :. j :. All)
            oldPJ j = trace ((show @String j) ++ "kiki : " ++ (show @String (computeS $ delay $ oldPJ1 j :: Array U DIM2 Double) ) ) $ oldPJ1 j
            -- oldPP i = getRows (knownIndices i) oldP --oldMu ! (Z :. j :. 0)
            xiP i = getRows (knownIndices i) (x i) -- !  (Z :. j :. 0)
            oldPP :: HasCallStack => Int -> Matrix D Double
            oldPP i = getRows (knownIndices i) oldP
            normOldPP i = sumAllS $ map (^2) $ oldPP i
            oldPPJ :: HasCallStack => Int -> Int -> Matrix D Double
            oldPPJ i j = extend (Any :. (1 :: Int)) $ (slice (oldPP i) (Any :. j) :: Vector D Double)
--            oldPPJ i j = trace ((show @String 89) ++ " : " ++ (show @String (computeS $ delay $ (oldPPJ1 i j) :: Array U DIM2 Double) ) ) $ oldPPJ1 i j
            w :: HasCallStack => Int -> Matrix D Double
            w i = map (+ (normOldPP i)) $ fromFunction (Z:. d3 :. d3) (\(Z :. y :. x) -> if y /= x then 0.0 else oldVariance) --mapDiagonal (+ oldVariance) $ (oldPP i) `mulS` (transpose $ oldPP i)
            invW i = trace ("2" :: String) $ invS $ w i
            sumPjXsubMu :: (HasCallStack) => Int -> Matrix D Double
            sumPjXsubMu i = sumListMatrices $ U.map (\j ->
              trace (("oldPPJ" :: String) ++ " : " ++ (show @String (computeS $ delay $ (oldPPJ i j) :: Array U DIM2 Double) ) ) $
                map (*((learnMatrix ! (Z:.j:.i)) - (oldMu ! (Z:.j:.0)))) (oldPJ j)) (knownIndices i) -- sumListMatrices $ U.map (\j -> map (* ((xij i j) - (oldMu j))) (oldP j)) knownIndices

            expX = fromFunction (Z :. d1 :. n)
                    (\(Z :. j :. i) -> if isNaN $ learnMatrix ! (Z :. j :. i) then (((oldPJ j) `mulS` (delay $ expT i)) ! (Z :. 0 :. 0)) + (oldMu ! (Z :. j :. 0)) else learnMatrix ! (Z :. j :. i))
            expXi i = extend (Any :. (1 :: Int)) $ slice expX (Any :. i)
            expXij i j = (expXi i) ! (Z :. j :. 0) -- if isNaN $ learnMatrix ! (Z :. j :. i) then (((oldPJ j) `mulS` (delay $ expT i)) ! (Z :. 0 :. 0)) + (oldMu ! (Z :. j :. 0)) else learnMatrix ! (Z :. j :. i)

            tiTrTi i =  map (* oldVariance) (invW i) +^ ((delay $ expT i) `mulS` (transpose $ expT i))

            expXiTrXi i = fromFunction (Z :. d3 :. d3)
              (\(Z :. k :. j) ->
                if | (isNaN $ learnMatrix ! (Z :. j :. i)) && (isNaN $ learnMatrix ! (Z :. k :. i)) && j /= k ->
                       (((delay $ (oldPJ j) `mulS` (expT i)) `mulS` (transpose $ oldPJ k)) ! (Z :. 0 :. 0))*oldVariance + (expXij i j)*(expXij i k)
                   | (isNaN $ learnMatrix ! (Z :. j :. i)) && (isNaN $ learnMatrix ! (Z :. k :. i)) && j == k ->
                        (((delay $ (oldPJ j) `mulS` (expT i)) `mulS` (transpose $ oldPJ k)) ! (Z :. 0 :. 0) + 1.0)*oldVariance + (expXij i j)*(expXij i k)
                   | (isNaN $ learnMatrix ! (Z :. j :. i)) && (not $ isNaN $ learnMatrix ! (Z :. k :. i)) ->
                        (expXij i j)*(learnMatrix ! (Z :. k :. i))
                   | (not $ isNaN $ learnMatrix ! (Z :. j :. i)) && (isNaN $ learnMatrix ! (Z :. k :. i)) ->
                        (learnMatrix ! (Z :. j :. i))*(expXij i k)
                   | (not $ isNaN $ learnMatrix ! (Z :. j :. i)) && (not $ isNaN $ learnMatrix ! (Z :. k :. i)) ->
                        (learnMatrix ! (Z :. j :. i))*(learnMatrix ! (Z :. k :. i)) )
            expXiTrTi i = fromFunction (Z :. d3 :. n)
              (\(Z :. k :. j) ->
                if   (isNaN $ learnMatrix ! (Z :. k :. j))
                then (((oldPJ k) `mulS` (delay $ invW i)) ! (Z :. k :. j))*(oldVariance) + ((expXi i) `mulS` (transpose $ expT i) ! (Z :. 0 :. 0))
                else ((x i) `mulS` (transpose $ expT i)) ! (Z :. 0 :. 0))

        in stepTwo ([x, expT, oldPJ, expXi, tiTrTi, expXiTrXi, expXiTrTi], knownIndices, oldP) iteration oldVariance

      stepTwo :: (HasCallStack) => ([Int -> Matrix D Double], Int -> [Int], Matrix D Double) -> Int -> Double -> (Matrix D Double, Double, Double)
      stepTwo ([x, expT, oldPJ, expXi, tiTrTi, expXiTrXi, expXiTrTi], knownIndices, oldP) iteration oldVariance =
        let newMu :: Matrix D Double
            newMu = trace ("lalala" ++ " : " ++ (show @String (computeS $ delay $ oldP :: Array U DIM2 Double) ) ) $
                      map (\x -> x/(fromIntegral n)) $ sumListMatrices $ U.map (\i -> (x i) -^ (oldP `mulS` (expT i))) [0..(n - 1)]
            newP = delay $ (sumListMatrices $ U.map (\i -> (expXiTrTi i) -^ (newMu `mulS` (transpose $ expT i))) [0..(n-1)])
                     `mulS` (delay $ invS $ sumListMatrices $ U.map (\i -> tiTrTi i) [0..(n-1)])
            varianceStep i = map (+ ((newMu `mulS` (transpose newMu) ! (Z :. 0 :. 0)) - 2*((newMu `mulS` (expXi i)) ! (Z :. 0 :. 0)) ) ) $
              (expXiTrXi i) -^ (map (*2) $ (expXiTrTi i) `mulS` (transpose $ newP))
                +^ (map (*(2 * (newMu `mulS` (transpose $ expT i)) ! (Z:.0:.0))) (transpose newP)) +^ (delay (newP `mulS` (tiTrTi i)) `mulS` (transpose $ newP))
            newVariance = (sum $ U.map (trace2S . computeS . varianceStep) [0..(n-1)])/(fromIntegral $ n*d3)
            diffVariance = trace ("23" :: String) $ abs $ newVariance - oldVariance
            maxDiffNewOldW = trace ("24" :: String) $ foldAllS max 0.0 $ map abs $ oldP -^ newP
            expLogLikelihood = 1.0
        in case stopParam of
              Left maxIteration -> if maxIteration > iteration then stepOne newP newVariance (iteration + 1) newMu else (newP,newVariance, expLogLikelihood)
              Right stopDiffBetweenIterations -> if (max diffVariance maxDiffNewOldW) > stopDiffBetweenIterations then stepOne newP newVariance (iteration + 1) newMu else (newP,newVariance, expLogLikelihood)
  in stepOne initMatrix initVariance 0 initMu
-}

convertPPCATypeSafeData
  :: (KnownNat y, KnownNat x)
  => (DimMatrix D y x Double, Double, Double)
  -> (Matrix D Double, Double, Double)
convertPPCATypeSafeData ((DimMatrix w), variance, expLogLikelihood) = (w, variance, expLogLikelihood)
