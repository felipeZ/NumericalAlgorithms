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


data WolfeData = WolfeData {
                            phi   ::   Double -> Double -- ^ Phi(ak)
                           ,phi'  ::  Double -> Double  -- ^ Phi derivative function
                           ,phi0  :: !Double            -- ^ Phi function evaluated at zero phi'(0)  
                           ,phi0' :: !Double            -- ^ Phi derivative evaluated at zero phi'(0)
                           ,wolfeP:: WolfeParameters    -- ^
                            }


-- | Implementation of the Line search Algorithm using the Strong Wolfe Condition as
-- | Described on: J. Nocedal, S. Wright, Numerical Optimization, second Edition, p. 60-62.
wolfeLineSearch :: Function
               -> FunGrad 
               -> Point    -- | Minimization direction 
               -> Point    -- | current point
               -> AlphaMax -- | Maximum Step length 
               -> Maybe WolfeParameters -- | These Parameters must fulfill 0 < C1 < C2 < 1
               -> Double 
wolfeLineSearch f gradF d xs aMax maybeWP = recWolfe wdata 0 a1 1                 
   where wdata     = WolfeData fun dervPhi' (f xs) (dot d $ gradF xs) wp  
         a1        = undefined
         wp        = fromMaybe (WP 1e-4 0.9) maybeWP
         fun       = f . evalArg
         funGrad   = gradF . evalArg 
         evalArg a = U.zipWith (+) xs $ U.map  (*a) d
         dervPhi'  = dot d . funGrad
         
recWolfe ::  WolfeData -> Double -> Double -> Int -> Double          
recWolfe wd@(WolfeData phi phi' phi0 phi0' (WP c1 c2) ) ak_1 ak step =
 if bool1 then action1 else if bool2 then action2 else if bool3 then action3 else action4
  
  where phi_k    = phi ak 
        c1akPhi' = c1 * ak * phi0'
        bool1    = (phi_k > phi0 + c1akPhi') ||              -- ak violates the sufficient decrease condition
                     (step > 1 && (phi ak >= phi ak_1 ))     -- or phi(ak) >= phi(ak_1)
        bool2    = (abs $ phi' ak) <= (negate $ c2 * phi0')  -- acceptable step length       
        bool3    = (phi' ak) >= 0                            -- the gradient points to the minimization direction
        action1  = zoom wd ak_1 ak
        action2  = ak
        action3  = zoom wd ak ak_1 
        action4  = recWolfe wd ak (choose ak) (succ step)
        choose   = undefined

-- | find the acceptable step length in the given interval
zoom :: WolfeData -> Double -> Double -> Double
zoom wd@(WolfeData phi phi' phi0 phi0' (WP c1 c2) ) alo ahi  =
   if bool1 then action1 else if bool2 then action2 else if bool3 then action3 else action4

  where aj = interpolate alo ahi
        phi_j    = phi aj 
        c1ajPhi' = c1 * aj * phi0'
        bool1    = (phi_j > phi0 + c1ajPhi') || (phi aj >= phi alo)
        bool2    = (abs $ phi' aj) <= (negate $ c2 * phi0')
        bool3    = ((phi' aj) *(ahi-alo)) >= 0 
        action1  = zoom wd alo aj 
        action2  = aj
        action3  = zoom wd aj alo
        action4  = zoom wd aj ahi

        
-- | Finds an alpha fulfilling the condtion phi(ak) <= phi(0) + c1 ak phi'(0) 
interpolate :: Double -> Double -> Double
interpolate = undefined
