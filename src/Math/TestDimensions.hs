{-# OPTIONS_GHC -fplugin GHC.TypeLits.KnownNat.Solver #-}
{-# OPTIONS_GHC -fplugin GHC.TypeLits.Normalise #-}

{-# LANGUAGE AllowAmbiguousTypes, FlexibleContexts, FlexibleInstances,
             MagicHash, MultiParamTypeClasses, PolyKinds,
             UndecidableInstances #-}

module Math.TestDimensions where

import Data.Array.Repa
import Data.Type.Bool
import qualified Data.Vec.Pull as V
import GHC.TypeLits hiding (someNatVal)
import GPLVM.Types
import GPLVM.Util hiding (trace2S)
import Numeric.Dimensions hiding (Head, Length, Tail)
import Numeric.Type.Evidence
import Numeric.TypedList hiding (All, Length, length, map)
import qualified Numeric.TypedList as TL hiding (Head, Length, Tail)
import Prelude (isNaN, log)

import Universum hiding (All, Any, Vector, map, natVal, toList, transpose, (++))
import qualified Universum as U

import Data.Array.Repa.Algorithms.Matrix
import Data.Singletons.Prelude (EnumFromTo)
import Data.Singletons.Prelude.List (Length)
import Data.Vinyl.TypeLevel (AllConstrained)

import GHC.Exts (unsafeCoerce#)
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import System.Random
import Unsafe.Coerce


--TypeFamily MulMatr Dim

newtype DimMatrix r (y :: Nat) (x :: Nat) a
  = DimMatrix { getInternal :: Matrix r a}

data PPCA = PPCA
  {  _learningData        :: Matrix D Double
   , desiredDimentions    :: Int
   , stopParameter        :: Either Int Double
   , _variance            :: Double
   , _W                   :: Matrix D Double
   , _finalExpLikelihood  :: Double
   }

type EvidenceList (c :: k -> Constraint) (xs :: [k])
  = TypedList (Evidence' c) xs

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
             withEvidence (inferPPCAInputMatrices @d @y1 @y2 @x2 @x1) $ convertPPCATypeSafeData $ func @d ld initM initVariance stopParameter
  in PPCA{..}

inferPPCAInputMatrices
  :: forall d y1 y2 x2 x1.
  AllConstrained KnownNat [d, y1, y2, x2, x1]
  => Evidence (x2 ~ d, y1 ~ y2, 1 <= x1)
inferPPCAInputMatrices
  | natVal (Proxy :: Proxy y1) /= natVal (Proxy :: Proxy y2) =
    error $ toText $ "dimentions y1 and y2 should be equal, but y1 = " U.++ (show (natVal (Proxy :: Proxy y1))) U.++ " and y2 = " U.++ (show (natVal (Proxy :: Proxy y2)))
  | natVal (Proxy :: Proxy x2) /= natVal (Proxy :: Proxy d) =
    error $ toText $ "dimentions x2 and d should be equal, but x2 = " U.++ (show (natVal (Proxy :: Proxy x2))) U.++ " and d = " U.++ (show (natVal (Proxy :: Proxy d)))
  | natVal (Proxy :: Proxy x1) < 1 =
    error $ toText ("Input matrix should have at least one column" :: String)
  | otherwise = unsafeCoerce# (E @((y1 ~ y1), (y1 ~ y1), (x1 <= x1)))

mulM
  :: forall y1 x1 y2 x2.
  ( AllConstrained KnownNat [x1, x2, y1, y2]
  , x1 ~ y2
  )
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> DimMatrix D y1 x2 Double
mulM (DimMatrix m1) (DimMatrix m2) = DimMatrix $ delay $ m1 `mulS` m2

emptyM
  :: forall y x. (KnownNat x, KnownNat y, (0 <=? x || 0 <=? y) ~ 'True)
  => DimMatrix D y x Double
emptyM =
  let x = fromIntegral $ natVal (Proxy :: Proxy x)
      y = fromIntegral $ natVal (Proxy :: Proxy y)
  in DimMatrix (delay $ fromListUnboxed (Z :. y :. x) []) :: DimMatrix D y x Double

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

(^++^) :: forall y1 x1 y2 x2.
  ( AllConstrained KnownNat [x1, x2, y1, y2]
  , y1 ~ y2
  )
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> DimMatrix D y2 (x1 + x2) Double
(^++^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 ++ m2

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

type family LessMax (x :: Nat) (yys :: [Nat]) :: Constraint where
  LessMax x '[]     = ()
  LessMax x (y ': ys) = ((y <= x), LessMax x ys)

withListOfIndexes
  :: forall y r . KnownNat y =>
  [Int] -> (forall (lng :: Nat). (lng <= y) => [Int] -> r) -> r
withListOfIndexes indexes f =
  let lng = fromIntegral $ length indexes
  in if (fromIntegral $ (natVal (Proxy :: Proxy y))) < maximum indexes
     then error "index is out of range"
     else case someNatVal (fromIntegral lng) of
      (SomeNat (Proxy :: Proxy tlng)) ->
        let ev = if (natVal (Proxy :: Proxy y)) <= lng
                 then error "list of indexes is too long"
                 else unsafeCoerce# (E :: Evidence (y ~ y)) :: Evidence (tlng <= y)
        in withEvidence ev $ f @tlng indexes

--data EitherTest a b = LeftT a | RightT b

type family AllLess (iis :: [Nat]) (max :: Nat) :: Constraint where
  AllLess '[] max = ()
  AllLess (i ': is) max = (i <= max, AllLess is max)

class MapExp (is :: [Nat]) (x1 :: Nat) (x2 :: Nat) where
  mapExp
    :: Proxy is
    -> ((LessThen x1) -> DimMatrix D x2 1 Double)
    -> DimMatrix D x2 (Length is) Double

instance (KnownNat x1,
          KnownNat x2,
          KnownNat i,
          KnownNat (Length is),
          AllLess iis x1,
          iis ~ (i ': is)) => MapExp (i ': is) x1 x2 where
  mapExp (_ :: Proxy iis) f = f (Less $ Proxy @i) ^++^ mapExp (Proxy @is) f

instance (KnownNat x1, KnownNat x2) => MapExp '[] x1 x2 where
  mapExp _ _ = emptyM :: DimMatrix D x2 0 Double

getColumn
  :: forall n x y. (AllConstrained KnownNat[n, y, x], n <= x)
  => DimMatrix D y x Double
  -> DimMatrix D y 1 Double
getColumn (DimMatrix m) =
  let i = fromIntegral $ natVal (Proxy :: Proxy n) :: Int
  in DimMatrix $ extend (Any :. (1 :: Int)) $ slice m (Any :. i)


mapColumn
  :: forall n x y0 y1 x1. (AllConstrained KnownNat[n, y0, y1, x], n <= x)
  => DimMatrix D y0 x Double
  -> (DimMatrix D y0 1 Double -> DimMatrix D y1 x1 Double)
  -> DimMatrix D y1 x1 Double
mapColumn dm f = f $ getColumn @n dm

withDeletedRows
  :: forall toDel x1 y1. (AllConstrained KnownNat [x1,y1,toDel], toDel <= y1)
  => DimMatrix D y1 x1 Double
  -> [Int]
  -> DimMatrix D (y1 - toDel) x1 Double
withDeletedRows (DimMatrix m) indexes = DimMatrix $
  deleteRows indexes m

toListM
  :: DimMatrix D x1 y1 Double
  -> [Double]
toListM (DimMatrix m) = toList m

type family LessPairs (is :: [Nat]) (n :: Nat) :: Constraint where
  LessPairs '[] n = ()
  LessPairs (i ': is) n = ((i <= n), (LessPairs is n))

{- withDelNumber
  :: forall (i :: Nat) (maxX :: Nat) (maxY :: Nat) r. (KnownNat i, i <= maxX)
  => (forall toDelT i0. (KnownNat toDelT, toDelT <= maxY, KnownNat i0, i0 <= maxX) => r)
  -> r
withDelNumber f = case someNatVal (fromIntegral $ toDel @i) of
  (SomeNat (Proxy :: Proxy toDelT)) ->
    let ev = if natVal (Proxy :: Proxy toDelT) > natVal (Proxy :: Proxy maxY)
             then error "list of indexes is too long"
             else unsafeCoerce# (E :: Evidence (i ~ i)) :: Evidence (toDelT <= maxY)
        in withEvidence ev $ f @toDelT @i -}

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

data LessThen (max :: Nat) = forall i. (KnownNat i, i <= max) => Less (Proxy i)

makeLessThenNat :: forall (max :: Nat). (KnownNat max) => Int -> LessThen max
makeLessThenNat i | i < 0                             = error "Natural number cannot be negative."
makeLessThenNat i | (natVal (Proxy :: Proxy max)) < (fromIntegral i) = error "Number can't be more then max type variable"
makeLessThenNat i | otherwise =
  case someNatVal (fromIntegral i) of
    (SomeNat (Proxy :: Proxy iT)) ->
      let ev = unsafeCoerce# (E :: Evidence (max <= max)) :: Evidence (iT <= max)
      in withEvidence ev $ Less (Proxy @iT)

mapAll :: [a -> r] -> [a] -> [r]
mapAll _ [] = []
mapAll [] _ = []
mapAll (f:fs) (x:xs) = f x : mapAll fs xs

--type family InferConstraint (complexConstraint :: Constraint) where
--  InferConstraint (iis ~ (i ': is)) = (Length iis) ~


emStepsMissed
  :: forall d y1 x1 y2 x2.
  ( HasCallStack
  , AllConstrained KnownNat [x1, x2, y1, y2]
  , x2 ~ d
  , y1 ~ y2
  , 1 <= x1
  )
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> Double
  -> Either Int Double
  -> (DimMatrix D y2 x2 Double, Double, Double, Maybe (DimMatrix D y1 x1 Double))
emStepsMissed learnMatrix initMatrix initVariance stopParam =
  let initMu = DimMatrix $ delay $ fromListUnboxed (Z:.d:.1) $ replicate d 0.0 :: DimMatrix D y1 1 Double
      n  = natVal (Proxy :: Proxy x1)
      d3 = natVal (Proxy :: Proxy y2)
      d  = fromIntegral $ natVal (Proxy :: Proxy y1)

      stepOne :: (HasCallStack)
              => DimMatrix D y1 1 Double
              -> DimMatrix D y2 x2 Double
              -> Double
              -> Int
              -> (DimMatrix D y2 x2 Double, Double, Double, Maybe (DimMatrix D y1 x1 Double))
      stepOne oldMu (!oldW) oldVariance iteration =
        let yi :: forall (i :: Nat). (KnownNat i, i <= x1) => DimMatrix D y1 1 Double
            yi = getColumn @i learnMatrix
            yList :: forall (i :: Nat). (KnownNat i, i <= x1) => [Double]
            yList = toListM $ yi @i
            unknownIndices :: forall (i :: Nat). (KnownNat i, i <= x1) => [Int]
            unknownIndices =
              let zipped = zip [0..] (yList @i)
              in U.map fst $ filter (isNaN . snd) zipped

            yP :: forall i toDel. (KnownNat i, KnownNat toDel, i <= x1, toDel <= y1) =>  DimMatrix D (y1 - toDel) 1 Double
            yP = withListOfIndexes @y1 (unknownIndices @i) (withDeletedRows (yi @i))

            muP :: forall i toDel. (KnownNat i, KnownNat toDel, i <= x1, toDel <= y1) => DimMatrix D (y1 - toDel) 1 Double
            muP = withListOfIndexes @y1 (unknownIndices @i) (withDeletedRows oldMu)

            oldWP :: forall i toDel. (KnownNat i, KnownNat toDel, i <= x1, toDel <= y1) => DimMatrix D (y1 - toDel) x2 Double
            oldWP = withListOfIndexes @y1 (unknownIndices @i) (withDeletedRows oldW)

            mP :: forall i toDel. (KnownNat i, KnownNat toDel, i <= x1, toDel <= y1) => DimMatrix D x2 x2 Double
            mP = mapDiagonalM (+ oldVariance) $ (transposeM (oldWP @i @toDel)) `mulM` (oldWP @i @toDel)

            invMP :: forall i toDel. (KnownNat toDel, KnownNat i, i <= x1, toDel <= y1, KnownNat (y2 - toDel)) => DimMatrix D x2 x2 Double
            invMP = invSM $ mP @i @toDel

            expXi :: (LessThen y2) -> (LessThen x1) -> DimMatrix D x2 1 Double
            expXi (Less (Proxy :: Proxy toDel)) (Less (Proxy :: Proxy i)) = ((invMP @i @toDel) `mulM` (transposeM (oldWP @i @toDel))) `mulM` ((yP @i @toDel) -^^ (muP @i @toDel))
--        in stepTwo (yP, expXi, invMP, unknownIndices, oldW) oldVariance iteration

--      stepTwo :: (HasCallStack) => ([Int -> Matrix D Double], Int -> [Int], Matrix D Double) -> Double -> Int -> (Matrix D Double, Double, Double, Maybe (Matrix D Double))
--      stepTwo (yP, expXi, invMP, unknownIndices, oldW) oldVariance iteration =
--        let -- expX = foldl1 (\acc i -> acc U.++ i) $ U.map expXi [0..(n-1)]
--            expX = reifyList unknownIndices
            ev1 :: forall (max :: Nat) (iis :: [Nat]) (i :: Nat) (is :: [Nat]). (1 <= max, iis ~ (EnumFromTo 0 (max -1)), iis ~ (i ': is)) => Evidence ((Length (EnumFromTo 0 (max -1))) ~ max)
            ev1 = unsafeCoerce# (E @(max ~ max))

--            ev :: forall (is :: [Nat]). (AllLess is x1, is ~ (EnumFromTo 0 (x1-1))) => Evidence (((Length is) ~ x1),(MapExp is x1 d))
--            ev = case (unsafeCoerce# (E @((x1 ~ x1),(x1 ~ x1),(x1 ~ x1))) :: Evidence (((Length is) ~ x1),(MapExp is x1 d), is ~ (EnumFromTo 0 (x1-1))))
--                   of E -> E
--            expX = withEvidence (ev1 @x1) $ mapExp (Proxy :: Proxy (EnumFromTo 0 (x1 - 1))) (withDelNumber expXi) :: DimMatrix D x2 x1 Double
                        {-
            mapExp
              :: [LessThen x1]
              -> [DimMatrix D x2 1 Double] -- (forall (toDel :: Nat) (i :: Nat). (KnownNat toDel, KnownNat i, i <= x1, toDel <= y2) => DimMatrix D y1 1 Double) -> DimMatrix D y1 (TL.Length is) Double
            mapExp ixs = U.map (withDelNumber expXi) ixs
-}
            toDel :: forall (i :: Nat). (KnownNat i, i <= x1) => Int
            toDel = length $ unknownIndices @i

            withDelNumber :: forall r. (LessThen y2 -> LessThen x1 -> r) -> LessThen x1 -> r
            withDelNumber f l@(Less (Proxy :: Proxy i)) = case someNatVal (fromIntegral $ toDel @i) of
                   (SomeNat (Proxy :: Proxy toDelT)) ->
                     let ev = if natVal (Proxy :: Proxy toDelT) > natVal (Proxy :: Proxy y1)
                              then error "list of indexes is too long"
                              else unsafeCoerce# (E :: Evidence (y2 ~ y2)) :: Evidence (toDelT <= y2)
                         toDelT = withEvidence ev $ Less $ (Proxy :: Proxy toDelT)
                     in f toDelT l

        in undefined
  in undefined
            {- newMu = extend (Any :. (1 :: Int)) $ meanColumnWithNan $ (learnMatrix -^ (oldW `mulS` expX))

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
              in transpose $ foldl1 (\acc j -> acc U.++ j) $ U.map (\x -> delay $ newWj x) [0..(d-1)]

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
-}
convertPPCATypeSafeData
  :: (KnownNat y, KnownNat x)
  => (DimMatrix D y x Double, Double, Double)
  -> (Matrix D Double, Double, Double)
convertPPCATypeSafeData ((DimMatrix w), variance, expLogLikelihood) = (w, variance, expLogLikelihood)
