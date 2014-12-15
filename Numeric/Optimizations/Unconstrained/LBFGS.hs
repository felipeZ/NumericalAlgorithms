{-# LANGUAGE BangPatterns, ViewPatterns #-}

module Numeric.Optimizations.Unconstrained.LBFGS (
              lBFGS 
              ) where

import Control.Exception (assert)
import Data.Sequence (
                     (<|)                    
                    ,ViewL(EmptyL , (:<) )
                    ,ViewR( (:>) )
                    ,viewl
                    ,viewr
                     )
import qualified Data.Sequence as S
import qualified Data.Vector.Unboxed as U
import Data.Array.Repa as R
import Data.Array.Repa.Algorithms.Matrix (
                        mmultS
                       ,transpose2S
                       )
import Text.Printf (printf)

-- Internal Modules
import  Numeric.Optimizations.TypesOptimization (
                          Function
                         ,FunGrad
                         ,Gradient
                         ,Iterations
                         ,Matrix
                         ,MaxSteps
                         ,Point
                         ,Step
                         ,Tolerance
                          )
  
import  Numeric.Optimizations.Tools (
             dot
            ,identity
            ,normVec
            ,scalarMatrix
            ,unboxed2Mtx
             )

import  Numeric.Optimizations.Unconstrained.WolfeCondition (wolfeLineSearch)

-- ==============> Data Types <==================

type StoreInfo = S.Seq (Point,Point) -- | Previous step and Gradients store to calculate the direction


-- ========================> <===============

lBFGS :: Monad m => Function   -- | Objective funxtion f(X) where X = {x1,x2...xn} 
                 -> FunGrad    -- | Objective function Gradient Df(X) = {df/dx1,..df/dxn} 
                 -> Point      -- | Initial point 
                 -> Matrix     -- | Initial Hessian Matrix
                 -> Iterations -- | Number of previous previous iteraction to store 
                 -> Tolerance  -- | Numerical Tolerance 
                 -> MaxSteps   -- | Maximum allowed steps
                 -> m (Either String Point)
lBFGS f gradF point guessHo mIterations delta maxSteps = recLBFGS point (Just guessHo) S.empty 1
  where recLBFGS !xs maybeHs info step =
          if (step > maxSteps)
             then return . Left $ printf "Convergence criterion not met after %d steps\n" step
             else do       
              let norma = normVec $ gradF xs        
              if (norma < delta)
                 then return $ Right xs
                 else do
                      let ho       = chooseH0 maybeHs info 
                          d        = calculateDir info ho (gradF xs) 
                          alpha    = wolfeLineSearch f gradF d xs 2 Nothing
                          alphaDir = U.map (*alpha) d 
                          newXs    = U.zipWith (+) xs alphaDir
                          g        = U.zipWith (-) (gradF newXs) (gradF xs)
                          s        = U.zipWith (-) newXs xs
                          newInfo  = updateInfo info (s,g) mIterations step
                      recLBFGS newXs Nothing newInfo (succ step) 


-- If it is the first step we use the provided Hessian Matrix otherwise we guess a new one
--  View Patterns is used to pattern match the head of the sequence 
chooseH0 :: Maybe Matrix -> StoreInfo -> Matrix  
chooseH0 (Just ho) _ = ho
chooseH0 Nothing (viewl -> (s,y) :< _ ) = scalarMatrix gamma $ identity dim 
  where gamma = (dot s y) * (recip $ dot y y)
        dim   = U.length s 

-- | Nocedal and Wright described this algorithm as a two-loop recursion one. A possible
-- | Implementation in Haskell  consist in a forward traverse of the store data structure
-- | then a backward traverse while accumulating.
calculateDir :: StoreInfo -- | previous stored m steps
             -> Matrix    -- | guess Hessian matrix
             -> Point     -- | gradient on current point
             -> Point     -- | minimization direction                   
calculateDir (viewl -> EmptyL) mtx grad  =  R.toUnboxed . mmultS mtx . unboxed2Mtx $ U.map negate grad
calculateDir info@(viewl -> (siU,yiU) :< _ ) mtx grad  = U.map negate $
  secondLoop info mtx . firstLoop info $ grad  
              
 where sM    = transpose2S $ unboxed2Mtx $ siU    
       yMT   = unboxed2Mtx $ yiU

-- |
firstLoop :: StoreInfo -> Point -> (Point,Point)
firstLoop info q0 = S.foldrWithIndex go (acc0, q0) info 
 where go i (si,yi) (acc,qi) = let ai     = calcAlpha si yi qi
                                   newQ   = U.zipWith (-) qi $ U.map (*ai) yi
                                   newAcc = U.cons ai acc 
                               in (newAcc,newQ)
       calcAlpha si yi qi = (rho si yi) * (dot si qi)
       rho si yi          = recip $ dot si yi   -- ( yT s) ^-1
       acc0               = U.empty
{-# INLINE firstLoop #-}

-- | 
secondLoop ::  StoreInfo -> Matrix -> (Point,Point) -> Point
secondLoop info guessH (as,q) =  S.foldrWithIndex go r0 $ S.reverse info       
      
 where go i (si,yi) ri = let b    = calcBeta si yi ri
                             ai   = as U.! i
                         in U.zipWith (+) ri $ U.map (*(ai - b)) si 
       dim = U.length r0                       
       r0 =  R.toUnboxed $ mmultS guessH (unboxed2Mtx q) 
       rho si yi          = recip $ dot si yi   -- ( yT s) ^-1
       calcBeta si yi ri = (rho si yi) * (dot yi ri)
{-# INLINE secondLoop #-}       

-- If the number of step is bigger than m drop last tuple and append new info
updateInfo :: StoreInfo -> (Point,Point) -> Iterations -> Step -> StoreInfo
updateInfo st new mIterations step | step <= mIterations = new <| st
                                   | otherwise = case viewr st of
                                                      (xs :> _) -> new <| xs

