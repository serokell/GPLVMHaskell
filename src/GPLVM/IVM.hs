-- | based on Informative vector machine implemented as follows:
-- https://github.com/lawrennd/ivm/blob/master/matlab/ivmCreate.m
-- See also
-- https://papers.nips.cc/paper/2240-fast-sparse-gaussian-process-methods-the-informative-vector-machine.pdf


module GPLVM.IVM
       ( IVMInput (..)
       , IVMOutput (..)
       , runIVM
       ) where

import Universum hiding (Vector)

import Prelude (log)

import Control.Lens (makeLenses)
import Data.Array.Repa

import Data.List (insert, (\\))

import GPLVM.Types

import Numeric.LinearAlgebra.Repa hiding (Matrix, Vector)

import Statistics.Distribution
import Statistics.Distribution.Normal (standard)

data IVMInput = IVMInput
  { _inputLength      :: Int              -- the number of input points
  , _desiredSparsity  :: Int              -- the number of active points
  , _covarianceMatrix :: Matrix D Double  -- covariance matrix
  , _inputMatrix      :: Matrix D Double
  , _biasParameter    :: Double           -- bias parameter (required for Gaussian noise)
  }

data IVMOutput = IVMOutput
  { _firstMatrix  :: Matrix D Double
  , _secondMatrix :: Matrix D Double
  , _thirdMatrix  :: Matrix D Double
  , _activeSet    :: [Int]
  }

makeLenses ''IVMInput
makeLenses ''IVMOutput

runIVM
  :: IVMInput
  -> IVMOutput
runIVM ivm =
  case sparsity <= len of
    True -> undefined
    False -> error errorMsg
  where
    sparsity = ivm ^. desiredSparsity
    len = ivm ^. inputLength
    bias = ivm ^. biasParameter
    covariance = ivm ^. covarianceMatrix
    inputList = [1..len]

    emptyVec :: Vector D Double
    emptyVec = ADelayed (Z :. sparsity) (const 0)
    
    initialPi :: Matrix D Double
    initialPi = diagD emptyVec
    initialA = diagD . takeDiagD $ covariance
    initialH = emptyVec
    initialM = emptyVec

    mMatrix
      :: Int
      -> Matrix D Double
    mMatrix ind = undefined

    updateCovariance
      :: Int
      -> Matrix D Double
    updateCovariance ind =
      covariance -^ mMatrix ind

    zet -- required for alpha
      :: Int
      -> Vector D Double
      -> Vector D Double
      -> Matrix D Double
      -> Double
    zet ind vec1 vec2 mat =
      index vec1 (Z :. ind) * (index vec2 (Z :. ind) + bias) /
      sqrt (1 + index mat (Z :. ind :. ind))

    alpha
      :: Int
      -> Vector D Double
      -> Vector D Double
      -> Matrix D Double
      -> Double
    alpha ind vec1 vec2 mat =
      index vec1 (Z :. ind) {- * normal distr -} /
      cumulativeNorm (index vec2 (Z :. ind)) *
      sqrt (index mat (Z :. ind :. ind))

    mu
      :: Int
      -> Vector D Double
      -> Vector D Double
      -> Vector D Double
      -> Matrix D Double
      -> Double
    mu ind vec1 vec2 vec3 mat =
      let alph = alpha ind vec1 vec2 mat in
      let frac = index vec3 (Z :. ind) + bias / 1 + index mat (Z :. ind :. ind) in
      alph * (alph + frac)

    mValue
      :: Int
      -> Vector D Double
      -> Vector D Double
      -> Vector D Double
      -> Matrix D Double
      -> Double
    mValue ind vec1 vec2 vec3 mat =
      index vec3 (Z :. ind) + alpha ind vec1 vec2 mat / mu ind vec1 vec2 vec3 mat


    delta   -- differential entropy
      :: Int
      -> Vector D Double
      -> Double
    delta ind vec =
      - 0.5 * log (1 - (index initialA (Z :. ind :. ind)) * index vec (Z :. ind))

    cumulativeNorm = cumulative standard
    errorMsg =
      "desired sparsity should be <= to the length of the input set"
