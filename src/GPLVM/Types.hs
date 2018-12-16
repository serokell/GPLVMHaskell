module GPLVM.Types where

import Universum

import           Control.Lens (makeLenses)

import           Numeric.Matrix (Matrix, MatrixElement)
import qualified Numeric.Matrix as M

type Distribution a = Matrix a -> a

type Mean a = Matrix a -> a


-------

zeroMean :: Num a => Mean a
zeroMean = const 0

newtype LatentSpacePoints a = LatentSpacePoints { _unLatentSpacePoints :: Matrix a }
    deriving Eq

makeLenses ''LatentSpacePoints

toLatentSpacePoints
    :: forall a. MatrixElement a
    => Matrix a
    -> LatentSpacePoints a
toLatentSpacePoints = LatentSpacePoints . M.transpose

newtype ObservedData a = ObservedData { _unObservedData :: Matrix a }

makeLenses ''ObservedData

data GaussianProcess a = GaussianProcess {
      _pointsGP :: Matrix a
    , _outputGP :: ObservedData a
    , _kernelGP :: Matrix a -> Matrix a -> Matrix a
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
    case (M.isSquare kernelMatrix && (kernelMatrix == M.transpose kernelMatrix)) of
        False -> Nothing
        True -> return kernelMatrix

data GaussianProcessLatentVariable a = GaussianProcessLatentVariable {
      _GPLVMObserved :: ObservedData a
    , _kerGPLVM :: Matrix a -> Matrix a -> Matrix a
    , _distrGPLVM ::  Distribution a
    }

makeLenses ''GaussianProcessLatentVariable
