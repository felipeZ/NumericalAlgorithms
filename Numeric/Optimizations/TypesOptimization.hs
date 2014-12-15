
module Numeric.Optimizations.TypesOptimization where

import qualified Data.Vector.Unboxed as U
import Data.Array.Repa as R


type Matrix    = Array U DIM2 Double 

type Gradient  = Array U DIM2 Double

type RepaVector= Array U DIM1 Double 

type Point     = U.Vector Double

type Function  = Point -> Double

type FunGrad   = Point -> Point  

type Tolerance = Double 

type MaxSteps  = Int

type Iterations = Int

type Step       = Int

type VecUnbox   = U.Vector Double

-- | Data Types 

data WolfeParameters = WP Double Double 
