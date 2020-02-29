module GPLVM.Types
       ( Covariance
       , Distribution
       , EigenVectors
       , InputObservations (..)
       , LatentSpacePoints (..)
       , Matrix (..)
       , Mean (..)
       , ObservedData (..)
       , Vector (..)
       , unLatentSpacePoints
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

type Matrix r a = Array r DIM2 a

type Vector r a = Array r DIM1 a

type EigenVectors a = Matrix D a

-- Real matrix, represented as unboxed vector.
-------

newtype LatentSpacePoints a = LatentSpacePoints { _unLatentSpacePoints :: Matrix D a }
    deriving Eq

makeLenses ''LatentSpacePoints

newtype InputObservations a = InputObservations { _unInputObs :: Vector D a }
    deriving Eq

newtype ObservedData a = ObservedData { _unObservedData :: Matrix D a }
    deriving Eq

makeLenses ''ObservedData
makeLenses ''InputObservations
