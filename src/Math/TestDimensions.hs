{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise #-}

{-# LANGUAGE AllowAmbiguousTypes, EmptyCase, FlexibleContexts,
             FlexibleInstances, LambdaCase, MagicHash, MultiParamTypeClasses,
             PolyKinds, UndecidableInstances #-}

module Math.TestDimensions where

import Data.List ((\\))
import Data.Singletons.Decide ((:~:)(..), Decision(..), (%~))
import Data.Type.Natural hiding ((+))

import GHC.TypeLits (ErrorMessage(Text), TypeError)
--import Data.Array.Repa hiding (Z)
import qualified Data.Array.Repa as R
import Data.Type.Bool
import qualified Data.Vec.Pull as V
--import GHC.TypeLits hiding (someNatVal, (<=))
import GPLVM.Types
import GPLVM.Util
--import Numeric.Dimensions hiding (Head, Length, Nat, Tail, (-))
import Numeric.Type.Evidence
import Numeric.TypedList hiding (All, Length, length, map)
import qualified Numeric.TypedList as TL hiding (Head, Length, Tail)
import Prelude (isNaN, log)

import Universum hiding
  (All, Any, Nat, One, Vector, map, natVal, toList, transpose, (%~), (++))
import qualified Universum as U hiding (Nat, One)

import Data.Array.Repa.Algorithms.Matrix hiding (trace2S)
import qualified Data.Singletons as S
import qualified Data.Singletons.Prelude as SP (EnumFromTo)

import Data.Vinyl.TypeLevel (AllConstrained)
import GHC.TypeNats (natVal)

import GHC.Exts (unsafeCoerce#)
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector, diag)
import System.Random
import Unsafe.Coerce


--TypeFamily MulMatr Dim

newtype DimMatrix r (y :: Nat) (x :: Nat) a
  = DimMatrix { getInternal :: Matrix r a}

showM :: DimMatrix R.D y x Double -> String
showM (DimMatrix m) = show (R.computeS m :: R.Array R.U R.DIM2 Double)

data PPCA = PPCA
  {  _noMissedData        :: Bool
   , _learningData        :: Matrix R.D Double
   , desiredDimensions    :: Int
   , stopParameter        :: Either Int Double
   , _variance            :: Double
   , _W                   :: Matrix R.D Double
   , _finalExpLikelihood  :: Double
   , _restoredMatrix      :: Maybe (Matrix R.D Double)
   }


withMat :: Matrix R.D Double -> (forall (x :: Nat) (y :: Nat). (SingI y, SingI x) => DimMatrix R.D x y Double -> k) -> k
withMat m f =
    let (R.Z  R.:.  y  R.:.  x) = R.extent m
    in
    case toSing (intToNat y) of
      SomeSing (sy :: Sing m) -> withSingI sy $
        case toSing (intToNat x) of
          SomeSing (sx :: Sing n) -> withSingI sx $ f (DimMatrix @R.D @m @n m)

makePPCATypeSafe
  :: RandomGen gen
  => Matrix R.D Double   -- ^ learning vectors
  -> Int               -- ^ desired dimension
  -> Either Int Double -- ^ number of iterations or
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

mulM
  :: forall y1 x1 y2 x2.
  (x1 ~ y2)
  => DimMatrix R.D y1 x1 Double
  -> DimMatrix R.D y2 x2 Double
  -> DimMatrix R.D y1 x2 Double
mulM (DimMatrix m1) (DimMatrix m2) = DimMatrix $ R.delay  $ m1 `mulS` m2

emptyM
  :: forall y x. (SingI y, SingI x) => (x ~ Zero)
  => DimMatrix R.D y x Double
emptyM =
  let x = fromIntegral $ toNatural (sing :: Sing x)
      y = fromIntegral $ toNatural (sing :: Sing y)
  in DimMatrix (R.delay  $ R.fromListUnboxed (R.Z  R.:.  y  R.:.  x) []) :: DimMatrix R.D y x Double

transposeM
  :: DimMatrix R.D y x Double
  -> DimMatrix R.D x y Double
transposeM (DimMatrix m) = DimMatrix $ R.transpose m

mapMM :: (Double -> Double)
  -> DimMatrix R.D y x Double
  -> DimMatrix R.D y x Double
mapMM f (DimMatrix m) =  DimMatrix $ R.map f m

mapDiagonalM :: (Double -> Double)
  -> DimMatrix R.D y x Double
  -> DimMatrix R.D y x Double
mapDiagonalM f (DimMatrix m) = DimMatrix $ mapDiagonal f m

invSM
  :: (y ~ x)
  => DimMatrix R.D y x Double
  -> DimMatrix R.D y x Double
invSM (DimMatrix m) = DimMatrix $ R.delay  $ invS m

substractMeanM
  :: DimMatrix R.D y x Double
  -> DimMatrix R.D y x Double
substractMeanM (DimMatrix m) = DimMatrix $ substractMean m

(+^^) :: forall y1 x1 y2 x2.
  ( x1 ~ x2
  , y1 ~ y2
  )
  => DimMatrix R.D y1 x1 Double
  -> DimMatrix R.D y2 x2 Double
  -> DimMatrix R.D y2 x2 Double
(+^^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 R.+^ m2

(^++^) :: forall y1 x1 y2 x2.
  ( y1 ~ y2
  )
  => DimMatrix R.D y1 x1 Double
  -> DimMatrix R.D y2 x2 Double
  -> DimMatrix R.D y2 (x1 :+: x2) Double
(^++^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 R.++ m2

(-^^) :: forall y1 x1 y2 x2.
  ( x1 ~ x2
  , y1 ~ y2
  )
  => DimMatrix R.D y1 x1 Double
  -> DimMatrix R.D y2 x2 Double
  -> DimMatrix R.D y2 x2 Double
(-^^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 R.-^ m2

(*^^) :: forall y1 x1 y2 x2.
  ( x1 ~ x2
  , y1 ~ y2
  )
  => DimMatrix R.D y1 x1 Double
  -> DimMatrix R.D y2 x2 Double
  -> DimMatrix R.D y2 x2 Double
(*^^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 R.*^ m2

(^!^) :: forall y x. DimMatrix R.D y x Double
      -> (LessThen y, LessThen x)
      -> Double
(^!^) (DimMatrix m) ((Less (Proxy :: Proxy y1)), (Less (Proxy :: Proxy x1))) =
  let y1 = natToInt $ fromSing (Sing :: Sing y1)
      x1 = natToInt $ fromSing (Sing :: Sing x1)
  in m R.! (R.Z R.:. y1 R.:. x1)

cholM
  :: DimMatrix R.D y x Double
  -> DimMatrix R.D y x Double
cholM (DimMatrix m) = DimMatrix $ R.delay  $ chol $ trustSym $ R.computeS m

solveM :: forall x1 y1 x2 y2. (y1 ~ y2, x1 ~ y2)
       => DimMatrix R.D y1 x1 Double
       -> DimMatrix R.D y2 x2 Double
       -> DimMatrix R.D y1 x2 Double
solveM (DimMatrix m1) (DimMatrix m2) = DimMatrix $ R.delay $ m1 `solveS` m2

sumAllSM
  :: DimMatrix R.D y x Double
  -> Double
sumAllSM (DimMatrix m) = R.sumAllS m

foldAllSM
  :: (Double -> Double -> Double)
  -> Double
  -> DimMatrix R.D y x Double
  -> Double
foldAllSM f initValue (DimMatrix m) = R.foldAllS f initValue m

detSM
  :: DimMatrix R.D y x Double
  -> Double
detSM (DimMatrix m) = detS m

diagM
  :: forall x y. (y ~ x)
  => DimMatrix R.D y x Double
  -> DimMatrix R.D y One Double
diagM (DimMatrix m) = DimMatrix $ diag m

trace2SM
  :: DimMatrix R.D y x Double
  -> Double
trace2SM (DimMatrix m) = trace2S $ R.computeS m

type family LessMax (x :: Nat) (yys :: [Nat]) :: Constraint where
  LessMax x '[]     = ()
  LessMax x (y ': ys) = ((y <= x) ~ 'True, LessMax x ys)

withListOfIndexes
  :: forall y r. (SingI y)
  => [Int] -> (LessThen y -> [Int] -> r) -> r
withListOfIndexes [] f =
  let lng = 0
  in f (makeLessThenNat @y lng) []
withListOfIndexes indexes f =
  let lng = fromIntegral $ length indexes
  in if (fromIntegral $ toNatural (sing :: Sing y)) < maximum indexes
     then error "index is out of range"
     else f (makeLessThenNat @y lng) indexes

type family AllLess (iis :: [Nat]) (max :: Nat) :: Constraint where
  AllLess  '[] max = ()
  AllLess  (i ': is) max = ((i <= max) ~ 'True, AllLess  is max)

{-type family Length (iis :: [k]) :: Nat where
  Length '[] = Zero
  Length (i ': is) = One :+: Length is -}

{-class MapExp (is :: [Nat]) (x1 :: Nat) (x2 :: Nat) where
  mapExp
    :: Proxy is
    -> ((LessThen x1) -> DimMatrix R.D x2 One Double)
    -> DimMatrix R.D x2 (Length is) Double

instance (AllLess iis x1,
          SingI x2,
          (One <= (Length iis)) ~ 'True,
          iis ~ (i ': is)) => MapExp (i ': is) x1 x2 where
  mapExp (_ :: Proxy iis) f = f (Less $ Proxy @i) ^++^ mapExp (Proxy @is) f

instance (SingI x2) => MapExp '[] x1 x2 where
  mapExp _ _ = emptyM :: DimMatrix R.D x2 Zero Double
-}
getColumn
  :: forall n x y. (((n <= x) ~ 'True), SingI n)
  => DimMatrix R.D y x Double
  -> DimMatrix R.D y One Double
getColumn (DimMatrix m) =
  let i = (fromIntegral $ toNatural (Sing :: Sing n)) :: Int
  in DimMatrix $ R.extend (R.Any   R.:.  (1 :: Int)) $ R.slice m (R.Any R.:.  i)

extendX :: forall n x y. (SingI n)
        => DimMatrix R.D y One Double
        -> DimMatrix R.D y n Double
extendX (DimMatrix m) =
  let n = natToInt $ fromSing (Sing :: Sing n) :: Int
  in DimMatrix $ (R.extend (R.Any R.:. n) $ R.slice m (R.Any   R.:.  (0 :: Int)))

getRow
  :: forall n x y. (((n <= y) ~ 'True), SingI n)
  => DimMatrix R.D y x Double
  -> DimMatrix R.D One x Double
getRow (DimMatrix m) =
  let j = (fromIntegral $ toNatural (Sing :: Sing n)) :: Int
  in DimMatrix $ R.extend (R.Any   R.:.  (1 :: Int)  R.:.  R.All ) $ R.slice m (R.Any   R.:.  (j :: Int)  R.:.  R.All )

mapColumn
  :: forall n x y0 y1 x1. (((n <= x) ~ 'True), SingI n)
  => DimMatrix R.D y0 x Double
  -> (DimMatrix R.D y0 One Double -> DimMatrix R.D y1 x1 Double)
  -> DimMatrix R.D y1 x1 Double
mapColumn dm f = f $ getColumn @n dm

meanColumnM
  :: forall y x. DimMatrix R.D y x Double
  -> DimMatrix R.D x One Double
meanColumnM (DimMatrix m) = DimMatrix $ R.extend (R.Any   R.:.  (1 :: Int)) $ meanColumn m

meanColumnWithNanM
  :: forall y x. DimMatrix R.D y x Double
  -> DimMatrix R.D y One Double
meanColumnWithNanM (DimMatrix m) = DimMatrix $ R.extend (R.Any   R.:.  (1 :: Int)) $ meanColumnWithNan m

sumListMatricesM
  :: [DimMatrix R.D y x Double]
  -> DimMatrix R.D y x Double
sumListMatricesM [] = error "Empty matrices list"
sumListMatricesM mss = foldl1 (\ms m -> ms +^^ m) mss

withDeletedRows
  :: forall toDel y1 x1. DimMatrix R.D y1 x1 Double
  -> LessThen y1
  -> [Int]
  -> DimMatrix R.D (y1 - toDel) x1 Double
withDeletedRows (DimMatrix m) (Less (Proxy :: Proxy i)) indexes = DimMatrix $
  deleteRows indexes m

withdeletedColumns
  :: forall toDel y1 x1. DimMatrix R.D y1 x1 Double
  -> LessThen x1
  -> [Int]
  -> DimMatrix R.D y1 (x1 - toDel) Double
withdeletedColumns (DimMatrix m) (Less (Proxy :: Proxy i)) indexes = DimMatrix $
  deleteColumns indexes m

toListM
  :: DimMatrix R.D x1 y1 Double
  -> [Double]
toListM (DimMatrix m) = R.toList m

type family LessPairs (is :: [Nat]) (n :: Nat) :: Constraint where
  LessPairs '[] n = ()
  LessPairs (i ': is) n = ((i <= n) ~ 'True, (LessPairs is n))

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

type family (n :: Nat) :+: (m :: Nat) where
  (:+:) n ('Z) = n
  (:+:) n ('S l) = ('S n) :+: l

type family (n :: Nat) :-: (m :: Nat) :: Nat where
  (:-:) ('Z) ('S _) = TypeError ('Text "lal")
  (:-:) ('Z) ('Z)   = ('Z)
  (:-:) k ('Z)      = k
  (:-:) ('S k) ('S l) = k :-: l

data LessThen (max :: Nat) = forall (i :: Nat). (((i <= max) ~ 'True), SingI i) => Less (Proxy i)

makeLessThenNat :: forall (max :: Nat). (SingI max) => Int -> LessThen max
makeLessThenNat i =
  case toSing (intToNat i) of
    (SomeSing (tt:: Sing (uu :: Nat))) ->
          case (tt %< (Sing :: Sing max)) of
            (Proved LS) -> withSingI tt $ Less (Proxy @uu)
            (Disproved _) -> error "the first number is not less then the second"

isLess :: forall (a :: Nat) (b :: Nat).(SingI (a <= b)) => Sing a -> Sing b -> Decision (a <= b :~: 'True)
isLess a b = ((Sing :: Sing (a<=b)) %~ (Sing :: Sing ('True)))

data (a :: k) :<: (b :: k) where
  LS :: forall k (a :: k) (b :: k). ((a <= b) ~ 'True) => a :<: b

class LSDecide k where
  (%<) :: forall (a :: k) (b :: k). Sing a -> Sing b -> Decision (a :<: b)
  infix 4 %<

instance LSDecide Nat where
  (%<) a b | Just l <- lessNat a b = Proved l
           | otherwise = Disproved (\_ -> error "the first number is not less then the second")

lessNat :: forall (a :: Nat) (b :: Nat). Sing a -> Sing b -> Maybe (a :<: b)
lessNat a b | toNatural a <= toNatural b = Just (unsafeCoerce# (LS :: Zero :<: One))
            | otherwise = Nothing


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
      d3 = demote @y2
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
            yP iP _ = withListOfIndexes @y1 (unknownIndices (Less iP)) (withDeletedRows @toDel (yi (Proxy :: Proxy i)))

            muP :: forall i toDel. ((i <= x1) ~ 'True, (toDel <= y1) ~ 'True, SingI i) => DimMatrix R.D (y1 - toDel) One Double
            muP = withListOfIndexes @y1 (unknownIndices (Less $ (Proxy :: Proxy i))) (withDeletedRows @toDel oldMu)

            oldWP :: forall i toDel. ((i <= x1) ~ 'True, (toDel <= y1) ~ 'True, SingI i)
                  => Proxy i -> Proxy toDel -> DimMatrix R.D (y1 - toDel) x2 Double
            oldWP iP _ = withListOfIndexes @y1 (unknownIndices (Less iP)) (withDeletedRows @toDel oldW)

            mP :: (LessThen x1) -> (LessThen y2) -> DimMatrix R.D x2 x2 Double
            mP iL@(Less iP@(Proxy :: Proxy iT)) toDelL@(Less toDelP@(Proxy :: Proxy toDelT)) =
              mapDiagonalM (+ oldVariance) $ (transposeM (oldWP iP toDelP)) `mulM` (oldWP iP toDelP)

            invMP :: (LessThen x1) -> (LessThen y2) -> DimMatrix R.D x2 x2 Double
            invMP i toDel = invSM $ mP i toDel

            expXi :: (LessThen y2) -> (LessThen x1) -> DimMatrix R.D x2 One Double
            expXi toDelL@(Less toDelP@(Proxy :: Proxy toDel)) iL@(Less iP@(Proxy :: Proxy i)) = ((invMP iL toDelL) `mulM` (transposeM (oldWP iP toDelP))) `mulM` ((yP iP toDelP) -^^ (muP @i @toDel))

        in stepTwo (expXi, invMP, unknownIndices, oldW) oldVariance iteration

--      stepTwo :: (HasCallStack) => ([Int -> Matrix R.D Double], Int -> [Int], Matrix R.D Double) -> Double -> Int -> (Matrix R.D Double, Double, Double, Maybe (Matrix R.D Double))
      stepTwo (expXi, invMP, unknownIndices, oldW) oldVariance iteration =
        let yi :: forall i. (((i <= x1) ~ 'True), SingI i) => Proxy i -> DimMatrix R.D y1 One Double
            yi _ = withSomeSing (demote @i) $ \s -> withSingI s $ getColumn @i learnMatrix

            yP :: forall i toDel. ((i <= x1) ~ 'True, (toDel <= y1) ~ 'True, SingI i)
               => Proxy i -> Proxy toDel -> DimMatrix R.D (y1 - toDel) One Double
            yP iP _ = withListOfIndexes @y1 (unknownIndices (Less iP)) (withDeletedRows @toDel (yi (Proxy :: Proxy i)))

            expX_ ::forall (i :: Nat). ((i <= x1) ~ 'True ) => Sing i  -> DimMatrix R.D x2 i Double
            expX_ SZ = emptyM :: DimMatrix R.D x2 Zero Double
            expX_ (SS l) = case lemma1 l (Sing :: Sing x1) of LS -> withSingI l $ (expX_ l) ^++^ ((withDelNumber expXi) (Less (Proxy :: Proxy (i :-: One))))

            expX :: DimMatrix R.D x2 x1 Double
            expX = case lemma2 (Sing :: Sing x1) (Sing :: Sing x1) of LS -> expX_ (Sing :: Sing x1)

            lemma2 :: forall (n :: Nat) (x :: Nat). (n ~ x) => Sing n -> Sing x -> (x :<: x)
            lemma2 (SZ) (SZ) = LS
            lemma2 (SS l) (SS k) = case lemma2 l k of LS -> LS

            lemma1 :: forall (n :: Nat) (x :: Nat). (('S n <= x) ~ 'True) => Sing n -> Sing x -> (n :<: x)
            lemma1 SZ _ = LS
            lemma1 (SS l) (SS k) = case lemma1 l k of LS -> LS

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
                  expXPj jP _ = withListOfIndexes @x1 (unknownColumns (Less jP)) (withdeletedColumns @toDel expX)

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
                          j = natToInt $ fromSing $ (Sing :: Sing j)

                  yPj :: forall j toDel. (((j <= y1) ~ 'True), (toDel <= x1) ~ 'True, SingI j)
                      => Proxy j -> Proxy toDel -> DimMatrix R.D One (x1 - toDel) Double
                  yPj jP _ = withListOfIndexes @x1 (unknownColumns (Less jP)) (withdeletedColumns @toDel $ yj @j)

                  expXtrXj :: forall j toDel. (((j <= y1) ~ 'True), (toDel <= x1) ~ 'True, SingI j)
                           => Proxy j -> Proxy toDel -> DimMatrix R.D x2 One Double
                  expXtrXj jP toDelP = (expXPj jP toDelP) `mulM` (transposeM $ mapMM (\x -> x - newMu ^!^ (Less jP,Less $ (Proxy :: Proxy 'Z))) (yPj jP toDelP))
                    where j = natToInt $ fromSing $ (Sing :: Sing j)

                  newWj :: LessThen x1 -> LessThen y1 -> DimMatrix R.D x2 One Double
                  newWj  (Less toDelP) (Less jP) = (gj jP toDelP) `solveM` (expXtrXj jP toDelP)

                  newW_ :: forall (j :: Nat). ((j <= y2) ~ 'True ) => Sing j  -> DimMatrix R.D x2 j Double
                  newW_ SZ = emptyM :: DimMatrix R.D x2 Zero Double
                  newW_ (SS l) = case lemma1 l (Sing :: Sing y2) of LS -> withSingI l $ (newW_ l) ^++^ ((withDelColumns newWj) (Less (Proxy :: Proxy (j :-: One))))
                    where j = natToInt $ fromSing (SS l)

              in transposeM $ case lemma2 (Sing :: Sing y2) (Sing :: Sing y1) of LS -> newW_ (Sing :: Sing y2)

            newVariance =
                let newWPi :: forall toDel. ((toDel <= y1) ~ 'True)
                          => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) x2 Double
                    newWPi iL _ = withListOfIndexes @y1 (unknownIndices iL) (withDeletedRows @toDel newW)

                    newMuP :: forall toDel. ((toDel <= y1) ~ 'True)
                           => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) One Double
                    newMuP iL _ = withListOfIndexes @y1 (unknownIndices iL) (withDeletedRows @toDel newMu)

                    varianceStep :: LessThen x1 -> LessThen y1 -> Double
                    varianceStep iL@(Less iP) toDelL@(Less toDelP) = sumAllSM $
                      (mapMM (^(2 :: Int)) ((yP iP toDelP) -^^ ((newWPi iL toDelP) `mulM` (expXi toDelL iL)) -^^ (newMuP iL toDelP)))
                      +^^ (mapMM (*oldVariance) $ diagM $ ((newWPi iL toDelP) `mulM` (invMP iL toDelL)) `mulM` (transposeM $ (newWPi iL toDelP)))

                    knownVars = foldAllSM (\acc x -> if isNaN x then acc else acc + 1.0) 0.0 learnMatrix

                    newVariance_ :: forall (i :: Nat). ((i <= x1) ~ 'True ) => Sing i -> Double
                    newVariance_ SZ = 0.0
                    newVariance_ (SS l) =
                      case lemma1 l (Sing :: Sing x1)
                      of LS -> withSingI l $ (newVariance_ l) + ((withDelNumber (flip varianceStep)) (Less (Proxy :: Proxy (i :-: One))))
                in (case lemma2 (Sing :: Sing x1) (Sing :: Sing x1) of LS -> newVariance_ (Sing :: Sing x1))/knownVars

            expLogLikelihood :: Double
            expLogLikelihood =
              let newWPi :: forall toDel. ((toDel <= y1) ~ 'True)
                         => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) x2 Double
                  newWPi iL _ = withListOfIndexes @y1 (unknownIndices iL) (withDeletedRows @toDel newW)

                  newMuP :: forall toDel. ((toDel <= y1) ~ 'True)
                         => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) One Double
                  newMuP iL _ = withListOfIndexes @y1 (unknownIndices iL) (withDeletedRows @toDel newMu)

                  yCenteredP :: forall toDel. ((toDel <= y1) ~ 'True)
                             => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) One Double
                  yCenteredP iL@(Less iP) toDelP = (yP iP toDelP) -^^ (newMuP iL toDelP)

                  yCenteredProd :: forall toDel. ((toDel <= y1) ~ 'True)
                                => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) (y1 - toDel) Double
                  yCenteredProd iL@(Less iP) toDelP = (yCenteredP iL toDelP) `mulM` (transposeM $ yCenteredP iL toDelP)

                  invMY :: forall toDel. ((toDel <= y1) ~ 'True)
                        => LessThen x1 -> Proxy toDel -> DimMatrix R.D (y1 - toDel) (y1 - toDel) Double
                  invMY iL toDelP = mapDiagonalM (+newVariance) $ (newWPi iL toDelP) `mulM` (transposeM $ newWPi iL toDelP)

                  knownIndices :: LessThen x1 -> [Int]
                  knownIndices iL = [0..(d - 1)] \\ (unknownIndices iL)

                  expLogLikelihoodStep :: LessThen x1 -> LessThen y1 -> Double
                  expLogLikelihoodStep iL@(Less iP) toDelL@(Less toDelP) =
                    (fromIntegral $ length $ knownIndices iL)*(log $ 2*pi) + (log $ detSM $ invMY iL toDelP) +
                    (trace2SM $ (invMY iL toDelP) `solveM` (yCenteredProd iL toDelP))

                  expLogLikelihood_ :: forall (i :: Nat). ((i <= x1) ~ 'True ) => Sing i -> Double
                  expLogLikelihood_ SZ = 0.0
                  expLogLikelihood_ (SS l) =
                    case lemma1 l (Sing :: Sing x1)
                    of LS -> withSingI l $ (expLogLikelihood_ l) + ((withDelNumber (flip expLogLikelihoodStep)) (Less (Proxy :: Proxy (i :-: One))))

              in (-1.0)*(case lemma2 (Sing :: Sing x1) (Sing :: Sing x1) of LS -> expLogLikelihood_ (Sing :: Sing x1))/2.0

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
