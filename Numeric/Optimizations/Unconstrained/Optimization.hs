
module Numeric.Optimizations.Unconstrained.Optimization
       (
       minimization
       ) where

import Data.Vector.Unboxed as U
import Data.Array.Repa as R


-- Internal Modules
import Numeric.Optimizations.Unconstrained.BFGS (bfgs)
import Numeric.Optimizations.Unconstrained.LBFGS (lBFGS)

import Numeric.Utilities.Tools (identity)

import Numeric.NumericTypes (
                          Function
                         ,FunGrad
                         ,Gradient
                         ,Matrix
                         ,MaxSteps
                         ,Threshold
                         ,VecUnbox
                          )

-- ===========> <================

-- | Maximum displacement 
type MaxDisp  = Maybe Double 

-- | Algorithm Used to minmize the target function
data Algorithm =  BFGS    
                | LBFGS  
                | CG      
                | Steepest MaxDisp

-- |
minimization :: Monad m => Function -> FunGrad -> VecUnbox -> Algorithm -> m (Either String VecUnbox)
minimization fun gradF xs algorithm =
  let hess0 = identity $ U.length xs
  in case algorithm of
       BFGS      ->  bfgs fun gradF xs hess0   1e-8 20
       LBFGS     ->  lBFGS fun gradF xs hess0 4 1e-8 20
       otherwise -> return $ Left "Sorry that's methods is still pending for implementation"





