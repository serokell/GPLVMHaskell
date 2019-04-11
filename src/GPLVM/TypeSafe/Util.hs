{-# LANGUAGE AllowAmbiguousTypes, FlexibleContexts, MagicHash #-}

-- | Utility functions for type-safe PPCA
-- This module contains mostly redefined functions for DimMatrix wrapper
-- around the REPA matrix
module GPLVM.TypeSafe.Util
  ( cholM
  , deleteColumnsM
  , deleteRowsM
  , detSM
  , diagM
  , eigSHM
  , emptyM
  , extendX
  , extendY
  , foldAllSM
  , getColumn
  , getFirstColumns
  , getRow
  , invSM
  , makeLessThenNat
  , mapDiagonalM
  , mapMM
  , meanColumnM
  , meanColumnWithNanM
  , meanCovSM
  , mulM
  , pinvSM
  , solveM
  , substractMeanM
  , sumAllSM
  , sumListMatricesM
  , toListM
  , trace2SM
  , transposeM
  , withMat
  , withVec
  , vecToMat
  , (+^^)
  , (^++^)
  , (-^^)
  , (*^^)
  , (^!^)
  , Decisions
  ) where

import Universum hiding (All, Any, Nat, One, Vector,
                         map, toList, transpose,
                         (%~), (++))

import GPLVM.Types
import GPLVM.TypeSafe.Types
import GPLVM.Util

import Data.Array.Repa
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector, diag)

import Data.Singletons.Decide (Decision(..))
import Data.Type.Natural hiding (Z)

-- | Creates new Matrix with dimensions in its type and passes it to continuation.
-- Also creates SinI instances for dimensions
withMat
  :: Matrix D Double
  -> (forall (x :: Nat) (y :: Nat). (SingI y, SingI x) => DimMatrix D x y Double -> k)
  -> k
withMat m f =
    let (Z  :.  y  :.  x) = extent m
    in
    case toSing (intToNat y) of
      SomeSing (sy :: Sing m) -> withSingI sy $
        case toSing (intToNat x) of
          SomeSing (sx :: Sing n) -> withSingI sx $ f (DimMatrix @D @m @n m)

-- | Vector cousin of withMat
withVec
  :: Vector D a
  -> (forall (n :: Nat). SingI n => DimVector D n a -> k)
  -> k
withVec vec f =
  case toSing (intToNat n) of
    SomeSing (sy :: Sing m) -> withSingI sy $ f (DimVector @D @m vec)
  where
    (Z :. n) = extent vec

vecToMat
  :: forall (n :: Nat) a.
  SingI n
  => DimVector D n a
  -> DimMatrix D n One a
vecToMat (DimVector vec) =
  DimMatrix $
  fromFunction (Z :. dim :. 1) (\(Z :. a :. b) -> index vec (Z :. a))
  where
    dim = natToInt $ demote @n :: Int

-- | Creates type-level number which is less then given Nat type
makeLessThenNat :: forall (max :: Nat). (SingI max) => Int -> LessThen max
makeLessThenNat i =
  case toSing (intToNat i) of
    (SomeSing (tt:: Sing (uu :: Nat))) ->
          case (tt %< (Sing :: Sing max)) of
            (Proved LS) -> withSingI tt $ Less (Proxy @uu)
            (Disproved _) -> error "the first number is not less then given type-level number"

-- | Creates empty column with desired number of rows
emptyM
  :: forall y x. (SingI y, SingI x) => (x ~ Zero)
  => DimMatrix D y x Double
emptyM =
  let x = fromIntegral $ toNatural (sing :: Sing x)
      y = fromIntegral $ toNatural (sing :: Sing y)
  in DimMatrix (delay  $ fromListUnboxed (Z  :.  y  :.  x) []) :: DimMatrix D y x Double

-- | Returns the column with given index

getColumn
  :: forall n x y. (((n <= x) ~ 'True), SingI n)
  => DimMatrix D y x Double
  -> DimMatrix D y One Double
getColumn (DimMatrix m) =
  let i = (fromIntegral $ toNatural (Sing :: Sing n)) :: Int
  in DimMatrix $ extend (Any   :.  (1 :: Int)) $ slice m (Any :.  i)

-- | Replicates the given column
extendX :: forall n y. (SingI n)
        => DimMatrix D y One Double
        -> DimMatrix D y n Double
extendX (DimMatrix m) =
  let n = natToInt $ fromSing (Sing :: Sing n) :: Int
  in DimMatrix $ (extend (Any :. n) $ slice m (Any   :.  (0 :: Int)))

-- | Replicates the given row
extendY :: forall n x. (SingI n)
        => DimMatrix D One x Double
        -> DimMatrix D n x Double
extendY (DimMatrix m) =
  let n = natToInt $ fromSing (Sing :: Sing n) :: Int
  in DimMatrix $ (extend (Any :. n :. All) $ slice m (Any   :.  (0 :: Int) :. All))

-- | Returns the row with given index
getRow
  :: forall n x y. (((n <= y) ~ 'True), SingI n)
  => DimMatrix D y x Double
  -> DimMatrix D One x Double
getRow (DimMatrix m) =
  let j = (fromIntegral $ toNatural (Sing :: Sing n)) :: Int
  in DimMatrix $ extend (Any   :.  (1 :: Int)  :.  All ) $ slice m (Any   :.  (j :: Int)  :.  All )

-- | Column of the mean values along every dimension
meanColumnM
  :: forall y x. DimMatrix D y x Double
  -> DimMatrix D x One Double
meanColumnM (DimMatrix m) = DimMatrix $ extend (Any   :.  (1 :: Int)) $ meanColumn m

-- | Similar to "meanColumnM" but it does not take into account NaN values
meanColumnWithNanM
  :: forall y x. DimMatrix D y x Double
  -> DimMatrix D y One Double
meanColumnWithNanM (DimMatrix m) = DimMatrix $ extend (Any   :.  (1 :: Int)) $ meanColumnWithNan m


-- todo: make sure that length of deleting indexes is equal to toDel
-- of change the list to list with length in its type
deleteRowsM
  :: forall toDel y1 x1. DimMatrix D y1 x1 Double
  -> [Int]
  -> DimMatrix D (y1 - toDel) x1 Double
deleteRowsM (DimMatrix m) indexes = DimMatrix $
  deleteRows indexes m

deleteColumnsM
  :: forall toDel y1 x1. DimMatrix D y1 x1 Double
  -> [Int]
  -> DimMatrix D y1 (x1 - toDel) Double
deleteColumnsM (DimMatrix m) indexes = DimMatrix $
  deleteColumns indexes m

getFirstColumns
  :: forall toSave y1 x1. (SingI toSave, toSave <= x1 ~ 'True)
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y1 toSave Double
getFirstColumns (DimMatrix m@(ADelayed (Z :. y :. x) f)) =
  let desiredColumns = natToInt $ fromSing (Sing :: Sing toSave)
  in DimMatrix $ ADelayed (Z :. y :. desiredColumns) f

sumListMatricesM
  :: [DimMatrix D y x Double]
  -> DimMatrix D y x Double
sumListMatricesM [] = error "Empty matrices list"
sumListMatricesM mss = foldl1 (\ms m -> ms +^^ m) mss

-- == Numeric functions from REPA of Numeric.LinearAlgebra for DimMatrix ==

mulM
  :: forall y1 x1 y2 x2.
  (x1 ~ y2)
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> DimMatrix D y1 x2 Double
mulM (DimMatrix m1) (DimMatrix m2) = DimMatrix $ delay  $ m1 `mulS` m2

transposeM
  :: DimMatrix D y x Double
  -> DimMatrix D x y Double
transposeM (DimMatrix m) = DimMatrix $ transpose m

mapMM :: (Double -> Double)
  -> DimMatrix D y x Double
  -> DimMatrix D y x Double
mapMM f (DimMatrix m) =  DimMatrix $ map f m

mapDiagonalM :: (Double -> Double)
  -> DimMatrix D y x Double
  -> DimMatrix D y x Double
mapDiagonalM f (DimMatrix m) = DimMatrix $ mapDiagonal f m

invSM
  :: (y ~ x)
  => DimMatrix D y x Double
  -> DimMatrix D y x Double
invSM (DimMatrix m) = DimMatrix $ delay  $ invS m

pinvSM
  :: DimMatrix D y x Double
  -> DimMatrix D x y Double
pinvSM (DimMatrix m) = DimMatrix $ delay $ pinvS m

substractMeanM
  :: DimMatrix D y x Double
  -> DimMatrix D y x Double
substractMeanM (DimMatrix m) = DimMatrix $ substractMean m

(+^^) :: forall y1 x1 y2 x2.
  ( x1 ~ x2
  , y1 ~ y2
  )
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> DimMatrix D y2 x2 Double
(+^^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 +^ m2

(^++^) :: forall y1 x1 y2 x2.
  ( y1 ~ y2
  )
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> DimMatrix D y2 (x1 :+: x2) Double
(^++^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 ++ m2

(-^^) :: forall y1 x1 y2 x2.
  ( x1 ~ x2
  , y1 ~ y2
  )
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> DimMatrix D y2 x2 Double
(-^^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 -^ m2

(*^^) :: forall y1 x1 y2 x2.
  ( x1 ~ x2
  , y1 ~ y2
  )
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> DimMatrix D y2 x2 Double
(*^^) (DimMatrix m1) (DimMatrix m2) = DimMatrix $ m1 *^ m2

(^!^) :: forall y x. DimMatrix D y x Double
      -> (LessThen y, LessThen x)
      -> Double
(^!^) (DimMatrix m) ((Less (Proxy :: Proxy y1)), (Less (Proxy :: Proxy x1))) =
  let y1 = natToInt $ fromSing (Sing :: Sing y1)
      x1 = natToInt $ fromSing (Sing :: Sing x1)
  in m ! (Z :. y1 :. x1)

cholM
  :: DimMatrix D y x Double
  -> DimMatrix D y x Double
cholM (DimMatrix m) = DimMatrix $ delay  $ chol $ trustSym $ computeS m

meanCovSM
  :: DimMatrix D y x Double
  -> (DimMatrix D One x Double, DimMatrix D x x Double)
meanCovSM (DimMatrix m) =
  let (meanVec, covMatrix) = meanCovS m
  in ((DimMatrix $ extend (Any  :. (1 :: Int) :. All) meanVec), DimMatrix $ delay $ unSym covMatrix)

solveM
  :: forall x1 y1 x2 y2.
  ( y1 ~ y2
  , x1 ~ y2
  )
  => DimMatrix D y1 x1 Double
  -> DimMatrix D y2 x2 Double
  -> DimMatrix D y1 x2 Double
solveM (DimMatrix m1) (DimMatrix m2) = DimMatrix $ delay $ m1 `solveS` m2

sumAllSM
  :: DimMatrix D y x Double
  -> Double
sumAllSM (DimMatrix m) = sumAllS m

foldAllSM
  :: (Double -> Double -> Double)
  -> Double
  -> DimMatrix D y x Double
  -> Double
foldAllSM f initValue (DimMatrix m) = foldAllS f initValue m

detSM
  :: DimMatrix D y x Double
  -> Double
detSM (DimMatrix m) = detS m

eigSHM
  :: (x ~ y, y1 ~ x) => DimMatrix D y x Double
  -> (DimMatrix D y1 One Double, DimMatrix D y x Double)
eigSHM (DimMatrix m) =
  let hermM = trustSymS m
      (eigenVals, eigenVecs) = eigSH hermM
  in ((DimMatrix $ extend (Any   :.  (1 :: Int)) eigenVals), DimMatrix $ delay eigenVecs)

diagM
  :: forall x y. (y ~ x)
  => DimMatrix D y x Double
  -> DimMatrix D y One Double
diagM (DimMatrix m) = DimMatrix $ diag m

trace2SM
  :: DimMatrix D y x Double
  -> Double
trace2SM (DimMatrix m) = trace2S $ computeS m

toListM
  :: DimMatrix D x1 y1 Double
  -> [Double]
toListM (DimMatrix m) = toList m

-- | this type family is introduced as a syntax-sugar for multiple decisions
type family Decisions k :: Type

type instance Decisions (a, b, c) = (Decision a, Decision b, Decision c)
