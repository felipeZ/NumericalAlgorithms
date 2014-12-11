{-# LANGUAGE BangPatterns #-}

module BFGS  (
              bfgs 
             ) where

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
            ,identity
            ,normVec
            ,scalarMatrix
            ,unboxed2Mtx
             )

import WolfeCondition (wolfeLineSearch)

-- ===============================> <========================

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
              if (norma < delta) then return $ Right xs
                                 else do
                                    let d        = mmultS hs . unboxed2Mtx $ U.map negate $ gradF xs
                                        xsM      = unboxed2Mtx xs
                                        alpha    = wolfeLineSearch f gradF (R.toUnboxed d) xs 2 Nothing
                                        alphaDir = scalarMatrix alpha d 
                                        newXs    = R.toUnboxed .computeUnboxedS $ xsM +^ alphaDir
                                        g        = U.zipWith (-) (gradF newXs) (gradF xs)
                                        s        = U.zipWith (-) newXs xs
                                    newHess <- updateHessian hs g s
                                    recBFGS newXs newHess (succ step) 

updateHessian :: Monad m =>
                    Matrix    -- | current Hessian Matrix
                 -> Point     -- | Delta Gradients GXj - GXj_1
                 -> Point     -- | Delta Points  Xj -Xj_1
                 -> m Matrix
updateHessian !mtx !gU !sU = computeUnboxedP $
                               hessk +^ rhoSST 
                                
 where g        = unboxed2Mtx gU
       s        = unboxed2Mtx sU
       rho      = recip $ dot gU sU
       ide      = identity $ U.length gU
       rhoSST   = scalarMatrix rho $ mmultS s sT
       rhoSYT   = scalarMatrix rho $ mmultS s gT
       rhoYST   = scalarMatrix rho $ mmultS g sT
       sT       = transpose2S s
       gT       = transpose2S g
       lft      = computeUnboxedS $ ide -^ rhoSYT
       rgt      = computeUnboxedS $ ide -^ rhoYST
       hessk    = mmultS lft $ mmultS mtx rgt 
{-# INLINE updateHessian #-}
