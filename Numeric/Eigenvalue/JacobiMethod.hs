{-# LANGUAGE BangPatterns #-}

module Numeric.Eigenvalue.JacobiMethod (
               jacobiP
              ) where

import Control.Arrow (second)
import Control.Applicative
import Control.Monad
import qualified Data.List as DL
import Data.Array.Repa         as R
import Data.Array.Repa.Unsafe  as R
import qualified Data.Vector.Unboxed as U


-- ==================> Internal Modules <===============
import Numeric.NumericTypes (Matrix, Step, Threshold, VecUnbox)
import Numeric.Utilities.Tools(identity, sortEigenData)

-- ============================> Types <=====================

data Parameters = Parameters !Double !Double !Double deriving Show

-- ============================> <======================
   
jacobiP :: Monad m => Matrix -> m (VecUnbox,Matrix)
jacobiP !arr = liftM (second (R.fromUnboxed sh) . sortEigenData .
                             second toUnboxed) $
                 loopJacobi arr ide 0 1.0e-9
  where sh@(Z :. dim :. _) = extent arr
        ide = identity dim

-- | Loop to carry out the corresponding rotation of the Jacobi Method
loopJacobi :: Monad m => Matrix -> Matrix -> Step -> Threshold -> m (VecUnbox,Matrix)
loopJacobi !arr !prr step tolerance
             | step > 5*dim^2         =  error "Jacobi method did not converge "
             | (arr ! mx) > tolerance = action
             | otherwise              = liftM2 (,) (toDiagonal arr) (return prr)
  where (Z:.dim:. _) = extent arr
        mx@(Z:.k:.l) = maxElemIndex arr
        aDiff        = (arr ! ix2 l l) - (arr ! ix2 k k)
        parameters   = calcParameters (arr ! mx) aDiff
        action       = do            
           newA <- rotateA arr mx parameters
           newP <- rotateP prr mx parameters
           loopJacobi newA newP (succ step) tolerance
           
                                         
rotateA :: Monad m => Matrix -> DIM2 -> Parameters -> m Matrix         
rotateA arr kl@(Z:. k :. l) (Parameters !s !t !tau) = computeUnboxedP $ unsafeTraverse arr id fun

  where fun f sh@(Z:.n:.m)  
                  | (n,m) == (k,l)           = 0.0
                  | (n,m) == (k,k)           = val - t*maxElem
                  | (n,m) == (l,l)           = val + t*maxElem
                  | n < k && m == k          = val - s*(gix2 n l + tau*val)
                  | n < k && m == l          = val + s*(gix2 n k - tau*val)
                  | k < m && m < l && n == k = val - s*(gix2 m l + tau*val)
                  | k < n && n < l && m == l = val + s*(gix2 k n - tau*val)
                  | m > l && n == k          = val - s*(gix2 l m + tau*val)
                  | m > l && n == l          = val + s*(gix2 k m - tau*val)
                  | otherwise = val 
 
               where gix2 a b = f $ ix2 a b
                     val      = f sh
                     maxElem  = f kl                     
{-# Inline rotateA #-}                     
                     
rotateP :: Monad m => Matrix -> DIM2 -> Parameters -> m Matrix         
rotateP prr kl@(Z:. k :. l) (Parameters !s !t !tau) = computeUnboxedP $ unsafeTraverse prr id  fun

  where fun f sh@(Z:.x:.y)                      
                 | y == k = val - s*(gix2 x l + tau*val)
                 | y == l = val + s*(gix2 x k - tau*val)
                 | otherwise = val
            where gix2 a b = f $ ix2 a b
                  val      = f sh
{-# Inline rotateP #-}                     

calcParameters ::  Double -> Double -> Parameters                 
calcParameters !maxElem !aDiff = Parameters s t tau
  where t = if abs maxElem < abs aDiff *1.0e-36 then maxElem/aDiff
            else let phi = aDiff/(2.0 *maxElem )
                     var = recip (abs phi + sqrt (1 + phi^2))
                 in if phi < 0 then negate var else var
        c = recip $ sqrt(t^2 + 1)
        s = t*c
        tau = s/(1 + c)
        
-- | Return the index of the largest off-diagonal element in the array
maxElemIndex  :: Matrix -> DIM2
maxElemIndex !arr  = R.fromIndex sh $ U.foldr fun 1 inds

  where inds = U.enumFromN 0 (dim^2) :: U.Vector Int
        sh@(Z:. dim :. _dim) = extent arr
        fun n acc= let sh2@(Z:. i:.j) = R.fromIndex sh n
                       sh3 = R.fromIndex sh acc
                   in if (i < j) && (abs (arr ! sh2) > abs (arr ! sh3)) then n else acc
       
toDiagonal :: Monad m => Matrix -> m VecUnbox        
toDiagonal mtx = liftM toUnboxed  $ computeUnboxedP $ unsafeBackpermute (ix1 dim) (\(Z:.i) -> ix2 i i) mtx
  where (Z:. dim :._) = extent mtx
