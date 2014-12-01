
module ConjugateGradients where

import Data.Vector.Unboxed as U
import Data.Array.Repa as R
import Data.Array.Repa.Algorithms.Matrix (
                       ,col
                       ,mmultP
                       ,mmultS
                       ,row
                       ,transpose2P
                       ,transpose2S
                       )

  
