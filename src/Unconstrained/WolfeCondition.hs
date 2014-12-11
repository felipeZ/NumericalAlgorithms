{-# LANGUAGE BangPatterns #-}

module WolfeCondition where

import Data.Maybe (fromMaybe)
import qualified Data.Vector.Unboxed as U

import Prelude as P

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
             )

-- Types
type AlphaMax  = Double


data WolfeData = WolfeData {
                            phi   ::  Double -> Double -- ^ Phi(ak)
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
wolfeLineSearch f gradF d xs aMax maybeWP = recWolfe wdata 0 a1 aMax 1                 
   where wdata     = WolfeData fun dervPhi' (f xs) (dot d $ gradF xs) wp  
         a1        = 1
         wp        = fromMaybe (WP 1e-4 0.9) maybeWP
         fun       = f . evalArg
         funGrad   = gradF . evalArg 
         evalArg a = U.zipWith (+) xs $ U.map  (*a) d
         dervPhi'  = dot d . funGrad
         
recWolfe ::  WolfeData -> Double -> Double -> Double -> Int -> Double
recWolfe wd@(WolfeData phi phi' phi0 phi0' (WP c1 c2) ) !ak_1 !ak aMax step =
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
        action4  = recWolfe wd ak choose aMax (succ step)
        choose   = interpolate wd ak aMax

-- | find the acceptable step length in the given interval
zoom :: WolfeData -> Double -> Double -> Double
zoom wd@(WolfeData phi phi' phi0 phi0' (WP c1 c2) ) alo ahi  =
   if bool1 then action1 else if bool2 then action2 else if bool3 then action3 else action4

  where aj       =   interpolate wd alo ahi
        phi_j    = phi aj 
        c1ajPhi' = c1 * aj * phi0'
        bool1    = (phi_j > phi0 + c1ajPhi') || (phi aj >= phi alo)
        bool2    = (abs $ phi' aj) <= (negate $ c2 * phi0')
        bool3    = ((phi' aj) *(ahi-alo)) >= 0 
        action1  = zoom wd alo aj 
        action2  = aj
        action3  = zoom wd aj alo        
        action4  =  zoom wd aj ahi

        
-- | Cubic interpolation of the guess Alpha_j in the interval
-- | bracketed by Alpha_Lower and Alpha_higher.
interpolate :: WolfeData -> Double -> Double -> Double
interpolate (WolfeData phi phi' phi0 phi0' (WP c1 c2) ) alo ahi  =
  if alo == 0 then ahi/2
              else (-b + sqr) / (3*a)
 where sqr = sqrt $ b^2 - 3*a * phi0'
       vs  = [(phi ahi) - phi0 - ahi*phi0',(phi alo) - phi0 - alo*phi0']
       mtx = [[alo^2, -ahi^2],[-alo^3,ahi^3]]
       fac = recip $ (alo^2)*(ahi^2)*(ahi-alo)
       [a,b] = P.map ((*fac) . P.sum . P.zipWith (*) vs) mtx
       
       
