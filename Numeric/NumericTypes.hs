
module Numeric.NumericTypes where 

-- ==============> Standard Modules and third party libraries <==========

import Data.Array.Repa          as R
import qualified Data.Vector.Unboxed as U

-- ==================> Internal Modules <===============

-- | Repa matrix representation
type Matrix = Array U DIM2 Double

type FunGrad   = Point -> Point  
type Function  = Point -> Double

type Gradient  = Array U DIM2 Double

-- | Repa matrix represnetation
type Matrix    = Array U DIM2 Double 

-- | Maximum number of cycles allowed
type MaxSteps  = Int

type RepaVector= Array U DIM1 Double 

-- | Iteration Step
type Step = Int

-- | Numerical tolerance
type Threshold = Double

-- | Unidimensional array of doubles representation
type VecUnbox  = U.Vector Double



-- =========> Data Types <============

data WolfeParameters = WP Double Double 
