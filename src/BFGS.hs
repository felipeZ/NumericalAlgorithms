
module BFGS where

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


-- Internal Modules
import TypesOptimization (
                          Gradient
                         ,Matrix
                         ,Vec
                          )
import Tools (
             dot
            ,scalarMatrix
             )

-- ===============================> <========================

-- initialGuess ::

updateHessian :: Monad m =>  Matrix -> Gradient -> Vec -> m Matrix
updateHessian mtx g s = computeUnboxedP $
                             mtx +^ gradTerm -^ hessTerm
                                
 where gradTerm = let gs = dot g s in scalarMatrix (1/gs) gradMtx 
       gradMtx  = mmultS g gT
       hessTerm = scalarMatrix (1/sTHs) $ mmultS hs hsT
       sTHs     = sumAllS $ mmultS sT hs 
       hs       = mmultS mtx s                
       sT       = transpose2S s
       gT       = transpose2S g
       hsT      = transpose2S hs
