module GPLVM.Types where

import Universum hiding (Vector)

import Control.Lens (makeLenses)

import Data.Array.Repa
import Data.Array.Repa.Repr.ForeignPtr
import Foreign.Storable
import GHC.TypeLits
import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)
import Numeric.Matrix (MatrixElement)
import qualified Numeric.Matrix as M

type Distribution r a = Matrix D r a -> a

type Mean r a = Matrix D r a -> a

type Covariance a = Herm a

newtype Matrix D r a = Matrix D (Array r DIM2 a)

newtype Vector r a = Vector (Array r DIM1 a)

type EigenVectors a = Matrix D F a

-- Real matrix, represented as unboxed vector.
-------

zeroMean :: Num a => Mean D a
zeroMean = const 0

newtype LatentSpacePoints a = LatentSpacePoints { _unLatentSpacePoints :: Matrix D D a }

makeLenses ''LatentSpacePoints

toLatentSpacePoints
    :: forall a. MatrixElement a
    => Matrix D a
    -> LatentSpacePoints a
toLatentSpacePoints = LatentSpacePoints . M.transpose

newtype ObservedData a = ObservedData { _unObservedData :: Matrix D a }

makeLenses ''ObservedData

data GaussianProcess a = GaussianProcess {
      _pointsGP :: Matrix D a
    , _outputGP :: ObservedData a
    , _kernelGP :: Matrix D a -> Matrix D a -> Matrix D a
    , _distrGP ::  Distribution a
    }

makeLenses ''GaussianProcess


covariaceMatrix
    :: forall a. MatrixElement a
    => GaussianProcess a
    -> Maybe (Matrix a)
covariaceMatrix gP = do
    let kernelFun = gP ^. kernelGP
    let points = gP ^. pointsGP
    let output = (gP ^. outputGP) ^. unObservedData
    let kernelMatrix = kernelFun points output
    case (M.isSquare kernelMatrix && (kernelMatrix == transpose kernelMatrix)) of
        False -> Nothing
        True -> return kernelMatrix

data GaussianProcessLatentVariable a = GaussianProcessLatentVariable {
      _GPLVMObserved :: ObservedData a
    , _kerGPLVM :: Matrix D a -> Matrix D a -> Matrix D a
    , _distrGPLVM ::  Distribution a
    }

makeLenses ''GaussianProcessLatentVariable -}
