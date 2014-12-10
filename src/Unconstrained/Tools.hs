
module Tools where

import qualified Data.Vector.Unboxed as U
import Data.Array.Repa as R
import Data.Array.Repa.Algorithms.Matrix (
                        col
                       ,mmultP
                       ,mmultS
                       ,row
                       ,transpose2P
                       ,transpose2S
                       )

-- Internal modules
import TypesOptimization (
                          Matrix
                         ,Point)

-- =============><==============

-- | scalar product between two vectors
dot :: Point -> Point -> Double
dot xs ys =  U.sum $ U.zipWith (*) xs ys

scalarMatrix :: Double -> Matrix -> Matrix
scalarMatrix s = computeUnboxedS . R.map (*s)

normVec :: Point -> Double
normVec xs = sqrt $ dot xs xs 

identity ::  Int -> Matrix
identity dim = computeUnboxedS $ fromFunction (ix2 dim dim) $
  \(Z:. x:. y) -> if x==y then 1 else 0

unboxed2Mtx :: U.Vector Double -> Matrix
unboxed2Mtx vs = R.fromUnboxed (ix2 dim 1) vs
  where dim = U.length vs
