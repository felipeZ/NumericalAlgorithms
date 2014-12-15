
module Numeric.Optimizations.Tools where

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
import Numeric.Optimizations.TypesOptimization (
                          Matrix
                         ,Point)

-- =============><==============

-- | scalar product between two vectors
dot :: Point -> Point -> Double
dot xs ys =  U.sum $ U.zipWith (*) xs ys

identity ::  Int -> Matrix
identity dim = computeUnboxedS $ fromFunction (ix2 dim dim) $
  \(Z:. x:. y) -> if x==y then 1 else 0

normVec :: Point -> Double
normVec xs = sqrt $ dot xs xs 

scalarMatrix :: Double -> Matrix -> Matrix
scalarMatrix s = computeUnboxedS . R.map (*s)

unboxed2Mtx :: U.Vector Double -> Matrix
unboxed2Mtx vs = R.fromUnboxed (ix2 dim 1) vs
  where dim = U.length vs

diagonal :: U.Vector Double -> Array U DIM2 Double
diagonal vs = computeUnboxedS $ traverse (identity dim) id
              $ \fun sh@(z :. i :. j) -> if i == j then vs U.! i
                                                   else fun sh
 where dim = U.length vs
                                                        
-- vector2Mtx :: U.Vector Double -> Matrix
-- vector2Mtx vs = R.fromUnboxed (ix2 1 dim ) vs
--   where dim = U.length vs
