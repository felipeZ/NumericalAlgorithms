
module Optimization where

import Data.Vector.Unboxed as U
import Data.Array.Repa as R


-- Internal Modules
import BFGS (bfgs)

import Tools (identity)

import TypesOptimization (
                          Function
                         ,FunGrad
                         ,Gradient
                         ,Matrix
                         ,MaxSteps
                         ,Point
                         ,Tolerance
                          )

-- ===========> <================

-- | Maximum displacement 
type MaxDisp  = Maybe Double 

-- | Algorithm Used to minmize the target function
data Algorithm =  BFGS    
                | L_BFGS  
                | CG      
                | Steepest MaxDisp

-- |
minimization :: Monad m => Function -> FunGrad -> Point -> Algorithm -> m (Either String Point)
minimization fun gradF xs algorithm =
  let hess0 = identity $ U.length xs
  in case algorithm of
       BFGS      -> bfgs fun gradF xs hess0 1e-6 100
       otherwise -> return $ Left "Sorry that's methods is still pending for implementation"



