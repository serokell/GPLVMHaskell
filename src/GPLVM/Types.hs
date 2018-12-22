module GPLVM.Types
       ( Covariance
       , Distribution
       , EigenVectors
       , InputObservations (..)
       , KernelFunction
       , LatentSpacePoints (..)
       , Matrix (..)
       , Mean (..)
       , ObservedData (..)
       , Vector (..)
       , unLatentSpacePoints
       , unVector
       , unMatrix
       , unObservedData
       , unInputObs
       ) where

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
    deriving Eq

makeLenses ''Matrix

newtype Vector r a = Vector { _unVector :: Array r DIM1 a}
    deriving Eq

makeLenses ''Vector

type EigenVectors a = Matrix D a

-- Real matrix, represented as unboxed vector.
-------

newtype LatentSpacePoints a = LatentSpacePoints { _unLatentSpacePoints :: Matrix D a }
    deriving Eq

makeLenses ''LatentSpacePoints


newtype InputObservations a = InputObservations { _unInputObs :: Matrix D a }
    deriving Eq

newtype ObservedData a = ObservedData { _unObservedData :: Matrix D a }
    deriving Eq

makeLenses ''ObservedData
makeLenses ''InputObservations

type KernelFunction a = Matrix D a -> Matrix D a -> Matrix D a
