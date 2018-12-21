module GPLVM.Types where

import Universum hiding (Vector, transpose)

import Control.Lens (makeLenses)

import Data.Array.Repa
import Data.Array.Repa.Repr.ForeignPtr
-- import Foreign.Storable
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

type Distribution r a = Matrix r a -> a

type Mean r a = Matrix r a -> a

type Covariance a = Herm a

newtype Matrix r a = Matrix { _unMatrix :: Array r DIM2 a }

makeLenses ''Matrix

newtype Vector r a = Vector (Array r DIM1 a)

type EigenVectors a = Matrix F a

-- Real matrix, represented as unboxed vector.
-------

zeroMean :: Num a => Mean D a
zeroMean = const 0

newtype LatentSpacePoints a = LatentSpacePoints { _unLatentSpacePoints :: Matrix D a }

makeLenses ''LatentSpacePoints

updateMatrix
    :: forall r a. (Array r DIM2 a -> Array r DIM2 a)
    -> Matrix r a
    -> Matrix r a
updateMatrix f (Matrix m) = Matrix (f m)

toLatentSpacePoints
    :: forall a. Matrix D a
    -> LatentSpacePoints a
toLatentSpacePoints = LatentSpacePoints . updateMatrix transpose

newtype ObservedData a = ObservedData { _unObservedData :: Matrix D a }

makeLenses ''ObservedData

data GaussianProcess a = GaussianProcess {
      _pointsGP :: Matrix D a
    , _outputGP :: ObservedData a
    , _kernelGP :: Matrix D a -> Matrix D a -> Matrix D a
    , _distrGP ::  Distribution D a
    }

makeLenses ''GaussianProcess


{-
covariaceMatrix
    :: forall a. GaussianProcess a
    -> Maybe (Matrix D a)
covariaceMatrix gP = do
    let kernelFun = gP ^. kernelGP
    let points = gP ^. pointsGP
    let output = (gP ^. outputGP) ^. unObservedData
    let kernelMatrix = kernelFun points output
    case (M.isSquare kernelMatrix && (kernelMatrix == transpose kernelMatrix)) of
        False -> Nothing
        True -> return kernelMatrix
-}

data GaussianProcessLatentVariable a = GaussianProcessLatentVariable {
      _GPLVMObserved :: ObservedData a
    , _kerGPLVM :: Matrix D a -> Matrix D a -> Matrix D a
    , _distrGPLVM ::  Distribution D a
    }

makeLenses ''GaussianProcessLatentVariable
