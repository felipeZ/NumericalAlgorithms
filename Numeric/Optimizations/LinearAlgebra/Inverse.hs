
module  Numeric.Optimizations.LinearAlgebra.Inverse (
               inverseP
                ) where


import Control.Monad 
import Data.Array.Repa         as R
import Data.Array.Repa.Unsafe  as R
import Data.Array.Repa.Algorithms.Matrix (
                        mmultP
                       ,mmultS
                       ,transpose2P
                       ,transpose2S
                       )
import qualified Data.Vector.Unboxed as U

-- ========> Internal Modules <===========

import Numeric.Optimizations.LinearAlgebra.Jacobi (jacobiP)
import Numeric.Optimizations.Tools  (diagonal)
import Numeric.Optimizations.TypesOptimization (Matrix)

-- =========> Types <============

-- =========> <============ 
inverseP :: Monad m => Matrix -> m Matrix
inverseP mtx = do
    (eigVals,vs) <- jacobiP mtx
    let diag  = diagonal $ U.map recip eigVals
        vsT   = transpose2S vs
    (mmultP vsT) <=< (mmultP diag) $ vs  
{-# INLINE inverseP #-}

        
