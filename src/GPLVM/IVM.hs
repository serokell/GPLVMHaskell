-- | based on Informative vector machine implemented as follows:
-- https://github.com/lawrennd/ivm/blob/master/matlab/ivmCreate.m
-- See also
-- https://papers.nips.cc/paper/2240-fast-sparse-gaussian-process-methods-the-informative-vector-machine.pdf


module GPLVM.IVM where

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
  { _inputLength      :: Int
  , _desiredSparsity  :: Int
  , _covarianceMatrix :: Matrix D Double
  , _inputMatrix      :: Matrix D Double
  , _biasParameter    :: Double
  }

data IVMOutput = IVMOutput
  { firstMatrix  :: Matrix D Double
  , secondMatrix :: Matrix D Double
  , thirdMatrix  :: Matrix D Double
  }

makeLenses ''IVMInput

runIVM
  :: IVMInput
  -> [Double]
runIVM ivm =
  case sparsity <= len of
    True -> do
      undefined
    False -> error errorMsg
  where
    sparsity = ivm ^. desiredSparsity
    len = ivm ^. inputLength
    bias = ivm ^. biasParameter
    covariance = ivm ^. covarianceMatrix
    inputList = [1..len]

    emptyVec = ADelayed (Z :. sparsity) (const 0)
    initialPi = diagD emptyVec
    initialA = diagD . takeDiagD $ covariance
    initialH = emptyVec
    initialM = emptyVec

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
      cumulativeNorm (index vec2 (Z :. ind)) * sqrt (index mat (Z :. ind :. ind))

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
