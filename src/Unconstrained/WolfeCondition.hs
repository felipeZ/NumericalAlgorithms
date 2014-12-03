{-# LANGUAGE BangPatterns #-}

module WolfeCondition where

import Data.Maybe (fromMaybe)
import qualified Data.Vector.Unboxed as U
import Data.Array.Repa as R
import Data.Array.Repa.Algorithms.Matrix (
                        mmultP
                       ,mmultS
                       ,transpose2P
                       ,transpose2S
                       )

-- Internal Modules
import TypesOptimization (
                          Function
                         ,FunGrad
                         ,Gradient
                         ,Matrix
                         ,MaxSteps
                         ,Point
                         ,Tolerance
                         ,WolfeParameters(..)
                          )
import Tools (
             dot
            ,normVec
            ,scalarMatrix
            ,unboxed2Mtx
             )

-- Types
type AlphaMax  = Double


-- | Implementation of the Line search Algorithm using the Strong Wolfe Condition as
-- | Described on: J. Nocedal, S. Wright, Numerical Optimization, second Edition, p. 60-62.
wolfeLineSearch :: Function
               -> FunGrad 
               -> Point    -- | Minimization direction 
               -> Point    -- | current point
               -> AlphaMax -- | Maximum Step length 
               -> Maybe WolfeParameters -- | These Parameters must fulfill 0 < C1 < C2 < 1
               -> Double 
wolfeLineSearch f gradF d xs aMax wp = recWolfe 0 a1 1  
               
   where a1        = undefined
         WP c1  c2 = fromMaybe (WP 1e-4 0.9) wp
         fun       = f . evalArg
         funGrad   = gradF . evalArg 
         evalArg a = U.zipWith (+) xs $ U.map  (*a) d
         phi0      = f xs               -- Phi function evaluated at zero phi'(0)  
         phi0'     = dot d $ gradF xs   -- Phi derivative evaluated at zero phi'(0)  
         phi' a    = dot d $ funGrad a   -- Phi derivative evaluated at ak  phi'(ak)  
         recWolfe  ak_1 ak step =
                   let phi_k    = fun ak
                       c1akPhi' = c1 * ak * phi0'
                       derv_Phi = dot d . gradF . evalArg $ ak
                       bool1    = (phi_k > phi0 + c1akPhi') || ((fun ak >= fun ak_1 ) &&  (step > 1))
                       bool2    = (abs $ phi' ak) <= (negate $ c2 * phi0')
                       bool3    = (phi' ak) >= 0
                       action1  = zoom ak_1 ak
                       action2  = ak
                       action3  = zoom ak ak_1
                   in undefined
         zoom      = undefined
