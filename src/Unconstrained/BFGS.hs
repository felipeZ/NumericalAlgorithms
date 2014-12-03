{-# LANGUAGE BangPatterns #-}

module BFGS where

import qualified Data.Vector.Unboxed as U
import Data.Array.Repa as R
import Data.Array.Repa.Algorithms.Matrix (
                        mmultP
                       ,mmultS
                       ,transpose2P
                       ,transpose2S
                       )
import Text.Printf (printf)

-- Internal Modules
import TypesOptimization (
                          Function
                         ,FunGrad
                         ,Gradient
                         ,Matrix
                         ,MaxSteps
                         ,Point
                         ,Tolerance
                          )
import Tools (
             dot
            ,normVec
            ,scalarMatrix
            ,unboxed2Mtx
             )

import WolfeCondition (wolfeLineSearch)

-- ===============================> <========================

-- initialGuess ::

bfgs :: Monad m => Function  -- | Objective funxtion f(X) where X = {x1,x2...xn} 
                -> FunGrad   -- | Objective function Gradient Df(X) = {df/dx1,..df/dxn} 
                -> Point     -- | Initial point 
                -> Matrix    -- | Initial Hessian Matrix
                -> Tolerance -- | Numerical Tolerance 
                -> MaxSteps  -- | Maximum allowed steps
                -> m (Either String Point)
bfgs f gradF point guessH delta maxSteps = recBFGS point guessH 0
  where recBFGS !xs !hs step =
          if (step > maxSteps)
             then return . Left $ printf "Convergence criterion not met after %d steps\n" step
             else do       
              let norma = normVec $ gradF xs        
              if (norma < delta) then return $ Right point
                                 else do
                                    let d        = mmultS hs . unboxed2Mtx $ U.map negate $ gradF xs
                                        xsM      = unboxed2Mtx xs
                                        alpha    = wolfeLineSearch f gradF (R.toUnboxed d) xs 1 Nothing
                                        alphaDir = scalarMatrix alpha d 
                                        newXs    = R.toUnboxed .computeUnboxedS $ xsM +^ alphaDir 
                                    newHess <- updateHessian hs (gradF xs) xs
                                    recBFGS newXs newHess (succ step) 

-- | Use the BFGS quasi-Newton method to update the minimization direction
updateHessian :: Monad m =>
                    Matrix    -- | current Hessian Matrix
                 -> Point     -- | Gradient on Point Xj
                 -> Point     -- | Point Xj
                 -> m Matrix
updateHessian !mtx !gU !sU = computeUnboxedP $
                               mtx +^ gradTerm -^ hessTerm
                                
 where g        = unboxed2Mtx gU
       s        = unboxed2Mtx sU
       gradTerm = let gs = dot gU sU in scalarMatrix (1/gs) gradMtx 
       gradMtx  = mmultS g gT
       hessTerm = scalarMatrix (1/sTHs) $ mmultS hs hsT
       sTHs     = sumAllS $ mmultS sT hs 
       hs       = mmultS mtx  s                
       sT       = transpose2S s
       gT       = transpose2S g
       hsT      = transpose2S hs
{-# INLINE updateHessian #-}


