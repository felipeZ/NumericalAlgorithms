
module Tools where

import Data.Vector.Unboxed as U
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
                         ,Vec)

-- =============><==============

-- | scalar product between two vectors
dot :: Vec -> Vec -> Double
dot xs ys =  sumAllS $ xs *^ ys

scalarMatrix :: Double -> Matrix -> Matrix
scalarMatrix s = computeUnboxedS . R.map (*s)
