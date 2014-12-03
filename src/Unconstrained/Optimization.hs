
module Optimization where

import Data.Vector.Unboxed as U
import Data.Array.Repa as R

  
-- ===========> <================

-- | Maximum displacement 
type MaxDisp  = Maybe Double 

-- | Algorithm Used to minmize the target function
data Algorithm =  BFGS    
                | L_BFGS  
                | CG      
                | Steepest MaxDisp

-- |
minimization :: Function -> FunGrad -> Point -> Algorithm -> Either String Point
minimzation
