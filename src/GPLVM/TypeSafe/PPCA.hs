{-# LANGUAGE AllowAmbiguousTypes, EmptyCase, FlexibleContexts, LambdaCase,
             MagicHash, MultiParamTypeClasses, PolyKinds,
             UndecidableInstances #-}

-- | In this module, two versions of PPCA are implemented using "type-safe"
-- matrix, which has its dimensions in type.
-- Due to that debugging of errors is much easier.
--
-- Implementation of the usual PPCA is easy
-- and the function looks very similar to the non-typesafe version
-- But if the changes in dimensions are depend from
-- runtime values(like in PPCA for missed values),
-- implementation is much more difficult and uses
-- singletons library with a noticeable performance overhead.

module GPLVM.TypeSafe.PPCA where

import GPLVM.Types
import GPLVM.TypeSafe.Types
import GPLVM.TypeSafe.Util
import GPLVM.Util

-- Universum
import Prelude (isNaN, log)
import Universum hiding
  (All, Any, Nat, One, Vector, map, toList, transpose, (%~), (++))
import qualified Universum as U

-- Common libraries
import Data.List ((\\))
import Data.Vinyl.TypeLevel (AllConstrained)
import System.Random

-- Type-level related libraries
import Data.Singletons.Decide ((:~:)(..), Decision(..), (%~))
import Data.Type.Natural

-- REPA
import qualified Data.Array.Repa as R

makePPCATypeSafe
  :: RandomGen gen
  => Matrix R.D Double   -- ^ learning vectors
  -> Int                 -- ^ desired dimension
  -> Either Int Double   -- ^ number of iterations or
                         --   the value of stop difference of the likelihood expectation between interations of EM.
  -> gen
  -> PPCA
makePPCATypeSafe leaningVectors desiredDimensions stopParameter generator =
  let _learningData@(R.ADelayed (R.Z  R.:.  n  R.:.  _) _) = leaningVectors
      initMatrix = randomMatrixD generator (n,desiredDimensions) :: Matrix R.D Double
      initVariance = fromIntegral $ fst . next $ generator :: Double
      _noMissedData = not $ isNaN $ runIdentity $ R.sumAllP _learningData
      (_W, _variance, _finalExpLikelihood, _restoredMatrix) =
        withMat _learningData $ \(ld :: DimMatrix R.D y1 x1 Double) ->
        withMat initMatrix $ \(initM :: DimMatrix R.D y2 x2 Double) ->
          case toSing (intToNat desiredDimensions) of
            (SomeSing (sDesired :: Sing d)) -> withSingI sDesired $
              let (ev1,ev2,ev3) = inferPPCAInputMatrices @d @y1 @y2 @x2 @x1
              in case ev1 of
                Disproved _ -> error $ toText @String $ "dimentions x2 and d1 should be equal, but they are not" -- y1 = " U.++ (show (natVal (Proxy :: Proxy y1))) U.++ " and y2 = " U.++ (show (natVal (Proxy :: Proxy y2)))
                Proved Refl -> case ev2 of
                  Disproved _ -> error $ toText @String $ "dimentions y1 and y2 should be equal, but they are not"
                  Proved Refl -> case ev3 of
                    Disproved _ -> error $ toText @String "Input matrix should have at least one column"
                    Proved LS -> if _noMissedData
                                 then convertPPCATypeSafeData $ emStepsFast @d ld initM initVariance stopParameter
                                 else convertPPCATypeSafeData $ emStepsMissed @d ld initM initVariance stopParameter
  in PPCA{..}

inferPPCAInputMatrices
  :: forall d y1 y2 x2 x1. (AllConstrained SingI '[d,y1,y2,x2,x1]) =>
     (Decision (x2 :~: d), Decision (y1 :~: y2), Decision (One :<: x1))
inferPPCAInputMatrices =
  let d1 = (sing :: Sing x2) %~ (sing :: Sing d)
      d2 = (sing :: Sing y1) %~ (sing :: Sing y2)
      d3 = (sing :: Sing One) %< (sing :: Sing x1)
  in (d1,d2,d3)


withListOfIndexes
  :: forall (y :: Nat) r. (SingI y)
  => [Int] -> ([Int] -> r) -> r
withListOfIndexes [] f = f []
withListOfIndexes indexes f =
  if (fromIntegral $ toNatural (sing :: Sing y)) < maximum indexes
  then error "index is out of range"
  else f indexes

emStepsFast
  :: forall d y1 x1 y2 x2.
  ( HasCallStack
  , x2 ~ d
  , y1 ~ y2
  , SingI x1
  , SingI y2
  )
  => DimMatrix R.D y1 x1 Double
  -> DimMatrix R.D y2 x2 Double
  -> Double
  -> Either Int Double
  -> (DimMatrix R.D y2 x2 Double, Double, Double, Maybe (DimMatrix R.D y1 x1 Double))
emStepsFast learnMatrix initMatrix initVariance stopParam =
  let learnMatrixCentered = transposeM $ substractMeanM (transposeM learnMatrix)
      n = natToInt $ demote @x1 :: Int
      d3 = natToInt $ demote @y2 :: Int

      stepOne :: (HasCallStack)
              => DimMatrix R.D y2 x2 Double
              -> Double
              -> Int
              -> (DimMatrix R.D y2 x2 Double, Double, Double, Maybe (DimMatrix R.D y1 x1 Double))
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
           => (DimMatrix R.D y12 x12 Double, DimMatrix R.D y22 x22 Double, DimMatrix R.D y32 x32 Double)
           -> DimMatrix R.D y42 x42 Double
           -> Double
           -> Int
           -> (DimMatrix R.D y52 x52 Double, Double, Double, Maybe (DimMatrix R.D y1 x1 Double))
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
              Left maxIteration -> if maxIteration > iteration then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood, Nothing)
              Right stopDiffBetweenIterations -> if (max diffVariance maxDiffNewOldW) > stopDiffBetweenIterations then stepOne newW newVariance (iteration + 1) else (newW,newVariance, expLogLikelihood, Nothing)

   in stepOne initMatrix initVariance 0

emStepsMissed
  :: forall d y1 x1 (y2 :: Nat) x2.
  ( HasCallStack
  , x2 ~ d
  , y1 ~ y2
  , (One <= x1) ~ 'True
  , SingI y2
  , SingI x1
  , SingI d
  )
  => DimMatrix R.D y1 x1 Double
  -> DimMatrix R.D y2 x2 Double
  -> Double
  -> Either Int Double
  -> (DimMatrix R.D y2 x2 Double, Double, Double, Maybe (DimMatrix R.D y1 x1 Double))
emStepsMissed learnMatrix initMatrix initVariance stopParam =
  let initMu = DimMatrix $ R.delay  $ R.fromListUnboxed (R.Z R.:. d R.:. 1) $ replicate d 0.0 :: DimMatrix R.D y1 One Double
      n  = demote @x1
      d  = natToInt $ demote @y1

      stepOne :: (HasCallStack)
              => DimMatrix R.D y1 One Double
              -> DimMatrix R.D y2 x2 Double
              -> Double
              -> Int
              -> (DimMatrix R.D y2 x2 Double, Double, Double, Maybe (DimMatrix R.D y1 x1 Double))
      stepOne oldMu (!oldW) oldVariance iteration =
        let yi :: forall i. (((i <= x1) ~ 'True), SingI i) => Proxy i -> DimMatrix R.D y1 One Double
            yi _ = withSomeSing (demote @i) $ \s -> withSingI s $ getColumn @i learnMatrix

            yList :: forall (i :: Nat). (((i <= x1) ~ 'True), SingI i) => [Double]
            yList = toListM $ yi (Proxy :: Proxy i)

            unknownIndices :: LessThen x1 -> [Int]
            unknownIndices (Less (Proxy :: Proxy i)) =
              let zipped = zip [0..] (yList @i)
              in U.map fst $ filter (isNaN . snd) zipped

            yP :: forall i toDel. ((i <= x1) ~ 'True, (toDel <= y1) ~ 'True, SingI i)
               => Proxy i -> Proxy toDel -> DimMatrix R.D (y1 - toDel) One Double
            yP iP _ = withListOfIndexes @y1 (unknownIndices (Less iP)) (deleteRowsM @toDel (yi (Proxy :: Proxy i)))

            muP :: forall i toDel. ((i <= x1) ~ 'True, (toDel <= y1) ~ 'True, SingI i) => DimMatrix R.D (y1 - toDel) One Double
            muP = withListOfIndexes @y1 (unknownIndices (Less $ (Proxy :: Proxy i))) (deleteRowsM @toDel oldMu)

            oldWP :: forall i toDel. ((i <= x1) ~ 'True, (toDel <= y1) ~ 'True, SingI i)
                  => Proxy i -> Proxy toDel -> DimMatrix R.D (y1 - toDel) x2 Double
            oldWP iP _ = withListOfIndexes @y1 (unknownIndices (Less iP)) (deleteRowsM @toDel oldW)

            mP :: (LessThen x1) -> (LessThen y2) -> DimMatrix R.D x2 x2 Double
            mP (Less iP@(Proxy :: Proxy iT)) (Less toDelP@(Proxy :: Proxy toDelT)) =
              mapDiagonalM (+ oldVariance) $ (transposeM (oldWP iP toDelP)) `mulM` (oldWP iP toDelP)

            invMP :: (LessThen x1) -> (LessThen y2) -> DimMatrix R.D x2 x2 Double
            invMP i toDel = invSM $ mP i toDel

            expXi :: (LessThen y2) -> (LessThen x1) -> DimMatrix R.D x2 One Double
            expXi toDelL@(Less toDelP@(Proxy :: Proxy toDel)) iL@(Less iP@(Proxy :: Proxy i)) = ((invMP iL toDelL) `mulM` (transposeM (oldWP iP toDelP))) `mulM` ((yP iP toDelP) -^^ (muP @i @toDel))

        in stepTwo (expXi, invMP, unknownIndices, oldW) oldVariance iteration


      stepTwo (expXi, invMP, unknownIndices, oldW) oldVariance iteration =
        let yi :: forall i. (((i <= x1) ~ 'True), SingI i) => Proxy i -> DimMatrix R.D y1 One Double
            yi _ = withSomeSing (demote @i) $ \s -> withSingI s $ getColumn @i learnMatrix

            yP :: forall i toDel. ((i <= x1) ~ 'True, (toDel <= y1) ~ 'True, SingI i)
               => Proxy i -> Proxy toDel -> DimMatrix R.D (y1 - toDel) One Double
            yP iP _ = withListOfIndexes @y1 (unknownIndices (Less iP)) (deleteRowsM @toDel (yi (Proxy :: Proxy i)))

            expX_ ::forall (i :: Nat). ((i <= x1) ~ 'True ) => Sing i  -> DimMatrix R.D x2 i Double
            expX_ SZ = emptyM :: DimMatrix R.D x2 Zero Double
            expX_ (SS l) = case lemma1 l (Sing :: Sing x1) of LS -> withSingI l $ (expX_ l) ^++^ ((withDelNumber expXi) (Less (Proxy :: Proxy (i - One))))

            expX :: DimMatrix R.D x2 x1 Double
            expX = case lemmaXLEqX (Sing :: Sing x1) (Sing :: Sing x1) of LS -> expX_ (Sing :: Sing x1)

            toDel :: forall (i :: Nat). (((i <= x1) ~ 'True), SingI i) => Proxy i -> Int
            toDel iP = length $ unknownIndices (Less iP)

            withDelNumber :: forall r. (SingI y2) => (LessThen y1 -> LessThen x1 -> r) -> LessThen x1 -> r
            withDelNumber f l@(Less iP) = f (makeLessThenNat @y2 (toDel iP)) l

            newMu :: DimMatrix R.D y2 One Double
            newMu = meanColumnWithNanM $ (learnMatrix -^^ (oldW `mulM` expX))

            newW :: DimMatrix R.D y2 x2 Double
            newW =
              let yj :: forall (j :: Nat). (((j <= y1) ~ 'True), SingI j) => DimMatrix R.D One x1 Double
                  yj = withSomeSing (demote @j) $ \s -> withSingI s $ getRow @j learnMatrix

                  yJList :: forall (j :: Nat). (((j <= y1) ~ 'True), SingI j) => [Double]
                  yJList = toListM $ yj @j

                  unknownColumns :: LessThen y1 -> [Int]
                  unknownColumns (Less (Proxy :: Proxy j)) =
                    let zipped = zip [0..] (yJList @j)
                    in U.map fst $ filter (isNaN . snd) zipped

                  knownColumns :: LessThen y1 -> [Int]
                  knownColumns j = [0..((natToInt n) - 1)] \\ (unknownColumns j)

                  expXPj :: forall j toDel. ((j <= y1) ~ 'True, (toDel <= x1) ~ 'True, SingI j)
                         => Proxy j -> Proxy toDel -> DimMatrix R.D x2 (x1 - toDel) Double
                  expXPj jP _ = withListOfIndexes @x1 (unknownColumns (Less jP)) (deleteColumnsM @toDel expX)

                  toDelColumns :: forall (j :: Nat). (((j <= y1) ~ 'True), SingI j) => Proxy j -> Int
                  toDelColumns jP = length $ unknownColumns (Less jP)

                  withDelColumns :: forall r. (SingI y2) => (LessThen x1 -> LessThen y1 -> r) -> LessThen y1 -> r
                  withDelColumns f l@(Less jP) = f (makeLessThenNat @x1 (toDelColumns jP)) l

                  sumInvMP :: forall j. ((j <= y1) ~ 'True, SingI j)
                           => Proxy j -> DimMatrix R.D x2 x2 Double
                  sumInvMP jP = sumListMatricesM $ U.map ((withDelNumber (flip invMP)) . (makeLessThenNat @x1)) (knownColumns $ Less $ jP)

                  gj :: forall j toDel. (((j <= y1) ~ 'True), (toDel <= x1) ~ 'True, SingI j)
                     => Proxy j -> Proxy toDel -> DimMatrix R.D x2 x2 Double
                  gj jP@(Proxy :: Proxy j) toDelP = (expXPJ `mulM` (transposeM expXPJ)) +^^ (mapMM (*oldVariance) (sumInvMP jP))
                    where expXPJ = expXPj jP toDelP

                  yPj :: forall j toDel. (((j <= y1) ~ 'True), (toDel <= x1) ~ 'True, SingI j)
                      => Proxy j -> Proxy toDel -> DimMatrix R.D One (x1 - toDel) Double
                  yPj jP _ = withListOfIndexes @x1 (unknownColumns (Less jP)) (deleteColumnsM @toDel $ yj @j)

                  expXtrXj :: forall j toDel. (((j <= y1) ~ 'True), (toDel <= x1) ~ 'True, SingI j)
                           => Proxy j -> Proxy toDel -> DimMatrix R.D x2 One Double
                  expXtrXj jP toDelP = (expXPj jP toDelP) `mulM` (transposeM $ mapMM (\x -> x - newMu ^!^ (Less jP,Less $ (Proxy :: Proxy 'Z))) (yPj jP toDelP))

                  newWj :: LessThen x1 -> LessThen y1 -> DimMatrix R.D x2 One Double
                  newWj  (Less toDelP) (Less jP) = (gj jP toDelP) `solveM` (expXtrXj jP toDelP)

                  newW_ :: forall (j :: Nat). ((j <= y2) ~ 'True ) => Sing j  -> DimMatrix R.D x2 j Double
                  newW_ SZ = emptyM :: DimMatrix R.D x2 Zero Double
                  newW_ (SS l) = case lemma1 l (Sing :: Sing y2) of LS -> withSingI l $ (newW_ l) ^++^ ((withDelColumns newWj) (Less (Proxy :: Proxy (j - One))))

              in transposeM $ case lemmaXLEqX (Sing :: Sing y2) (Sing :: Sing y1) of LS -> newW_ (Sing :: Sing y2)

            newVariance =
                let newWPi :: forall toDel. ((toDel <= y1) ~ 'True)
                          => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) x2 Double
                    newWPi iL _ = withListOfIndexes @y1 (unknownIndices iL) (deleteRowsM @toDel newW)

                    newMuP :: forall toDel. ((toDel <= y1) ~ 'True)
                           => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) One Double
                    newMuP iL _ = withListOfIndexes @y1 (unknownIndices iL) (deleteRowsM @toDel newMu)

                    varianceStep :: LessThen x1 -> LessThen y1 -> Double
                    varianceStep iL@(Less iP) toDelL@(Less toDelP) = sumAllSM $
                      (mapMM (^(2 :: Int)) ((yP iP toDelP) -^^ ((newWPi iL toDelP) `mulM` (expXi toDelL iL)) -^^ (newMuP iL toDelP)))
                      +^^ (mapMM (*oldVariance) $ diagM $ ((newWPi iL toDelP) `mulM` (invMP iL toDelL)) `mulM` (transposeM $ (newWPi iL toDelP)))

                    knownVars = foldAllSM (\acc x -> if isNaN x then acc else acc + 1.0) 0.0 learnMatrix

                    newVariance_ :: forall (i :: Nat). ((i <= x1) ~ 'True ) => Sing i -> Double
                    newVariance_ SZ = 0.0
                    newVariance_ (SS l) =
                      case lemma1 l (Sing :: Sing x1)
                      of LS -> withSingI l $ (newVariance_ l) + ((withDelNumber (flip varianceStep)) (Less (Proxy :: Proxy (i - One))))
                in (case lemmaXLEqX (Sing :: Sing x1) (Sing :: Sing x1) of LS -> newVariance_ (Sing :: Sing x1))/knownVars

            expLogLikelihood :: Double
            expLogLikelihood =
              let newWPi :: forall toDel. ((toDel <= y1) ~ 'True)
                         => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) x2 Double
                  newWPi iL _ = withListOfIndexes @y1 (unknownIndices iL) (deleteRowsM @toDel newW)

                  newMuP :: forall toDel. ((toDel <= y1) ~ 'True)
                         => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) One Double
                  newMuP iL _ = withListOfIndexes @y1 (unknownIndices iL) (deleteRowsM @toDel newMu)

                  yCenteredP :: forall toDel. ((toDel <= y1) ~ 'True)
                             => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) One Double
                  yCenteredP iL@(Less iP) toDelP = (yP iP toDelP) -^^ (newMuP iL toDelP)

                  yCenteredProd :: forall toDel. ((toDel <= y1) ~ 'True)
                                => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) (y1 - toDel) Double
                  yCenteredProd iL toDelP = (yCenteredP iL toDelP) `mulM` (transposeM $ yCenteredP iL toDelP)

                  invMY :: forall toDel. ((toDel <= y1) ~ 'True)
                        => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) (y1 - toDel) Double
                  invMY iL toDelP = mapDiagonalM (+newVariance) $ (newWPi iL toDelP) `mulM` (transposeM $ newWPi iL toDelP)

                  knownIndices :: LessThen x1 -> [Int]
                  knownIndices iL = [0..(d - 1)] \\ (unknownIndices iL)

                  expLogLikelihoodStep :: LessThen x1 -> LessThen y1 -> Double
                  expLogLikelihoodStep iL (Less toDelP) =
                    (fromIntegral $ length $ knownIndices iL)*(log $ 2*pi) + (log $ detSM $ invMY iL toDelP) +
                    (trace2SM $ (invMY iL toDelP) `solveM` (yCenteredProd iL toDelP))

                  expLogLikelihood_ :: forall (i :: Nat). ((i <= x1) ~ 'True ) => Sing i -> Double
                  expLogLikelihood_ SZ = 0.0
                  expLogLikelihood_ (SS l) =
                    case lemma1 l (Sing :: Sing x1)
                    of LS -> withSingI l $ (expLogLikelihood_ l) + ((withDelNumber (flip expLogLikelihoodStep)) (Less (Proxy :: Proxy (i - One))))

              in (-1.0)*(case lemmaXLEqX (Sing :: Sing x1) (Sing :: Sing x1) of LS -> expLogLikelihood_ (Sing :: Sing x1))/2.0

            restoredData =
              let muX :: DimMatrix R.D x2 One Double
                  muX = meanColumnM (transposeM expX)

                  finalMu = (newMu +^^ (newW `mulM` muX))

                  wtrw = (transposeM newW) `mulM` newW

                  factor1 = newW `mulM` (invSM wtrw)

                  factor2 = mapDiagonalM (+newVariance) wtrw

                  restoredCentered =  factor1 `mulM` factor2 `mulM` expX

              in restoredCentered +^^ extendX @x1 finalMu

            diffVariance = abs $ newVariance - oldVariance
            maxDiffNewOldW = foldAllSM max 0.0 $ mapMM abs $ oldW -^^ newW
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

convertPPCATypeSafeData
  :: (DimMatrix R.D y x Double, Double, Double, Maybe (DimMatrix R.D y1 x1 Double))
  -> (Matrix R.D Double, Double, Double, Maybe (Matrix R.D Double))
convertPPCATypeSafeData ((DimMatrix w), variance, expLogLikelihood, (Just (DimMatrix i))) = (w, variance, expLogLikelihood, Just i)
convertPPCATypeSafeData ((DimMatrix w), variance, expLogLikelihood, Nothing) = (w, variance, expLogLikelihood, Nothing)
