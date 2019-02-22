{-# LANGUAGE AllowAmbiguousTypes #-}

module Math.TestDimensions where

import Data.Array.Repa hiding ((++))
import GHC.TypeLits hiding (someNatVal)
import GPLVM.Types
import GPLVM.Util
import Numeric.Dimensions

import Prelude (log)

import Universum hiding (All, Any, Vector, map, natVal, toList, transpose)

import Data.Array.Repa.Algorithms.Matrix
import Data.Vector.Unboxed.Base (Unbox)
import Data.Vinyl.TypeLevel (AllConstrained)

import Foreign.Storable (Storable)

import Math.Matrix hiding (randomMatrixD)

import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import System.Random
import Unsafe.Coerce


--TypeFamily MulMatr Dim

data PPCA = PPCA
  {  _learningData       :: Matrix D Double
  , desiredDimentions    :: Int
  , stopParameter        :: Either Int Double
  , _variance            :: Double
  , _W                   :: Matrix D Double
  , _finalExpLikelihood  :: Double
  }

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
makePPCATypeSafe leaningVectors desiredDimensions stopParameter func generator =
  let _learningData@(ADelayed (Z :. n :. _) _) = leaningVectors
      initMatrix = randomMatrixD generator (n, desiredDimensions) :: Matrix D Double
      initVariance = fromIntegral $ fst . next $ generator :: Double
      (_W, _variance, _finalExpLikelihood) =
        withMat _learningData $ \(ld :: DimMatrix D y1 x1 Double) ->
        withMat initMatrix $ \(initM :: DimMatrix D y2 x2 Double) ->
        case someNatVal $ fromIntegral desiredDimensions of
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

convertPPCATypeSafeData
  :: (KnownNat y, KnownNat x)
  => (DimMatrix D y x Double, Double, Double)
  -> (Matrix D Double, Double, Double)
convertPPCATypeSafeData ((DimMatrix w), variance, expLogLikelihood) = (w, variance, expLogLikelihood)
