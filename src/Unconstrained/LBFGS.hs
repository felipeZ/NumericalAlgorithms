{-# LANGUAGE BangPatterns, ViewPatterns #-}

module LBFGS (
              lBFGS 
              ) where

import Data.Sequence (
                     (<|)
                    ,ViewL( (:<) )
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
import TypesOptimization (
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
import Tools (
             dot
            ,identity
            ,normVec
            ,scalarMatrix
            ,unboxed2Mtx
             )

import WolfeCondition (wolfeLineSearch)

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
chooseH0 _ (viewl -> (s,y) :< _ ) = scalarMatrix p $ mmultS sT yM
  where sT = transpose2S $ unboxed2Mtx s
        yM = unboxed2Mtx y
        p = recip $ dot y y                                     

-- | Nocedal and Wright described this algorithm as a two-loop recursion one. A possible
-- | Implementation in Haskell traverse consist in a forward traverse of the store data structure
-- | then a backward traverse while accumulating.
calculateDir :: StoreInfo -- | previous stored m steps
             -> Matrix    -- | guess Hessian matrix
             -> Point     -- | gradient on current point
             -> Point     -- | minimization direction                   
calculateDir = undefined
--foldrWithIndex :: (Int -> a -> b -> b) -> b -> Seq a -> b


-- If the number of step is bigger than m drop last tuple and append new info
updateInfo :: StoreInfo -> (Point,Point) -> Iterations -> Step -> StoreInfo
updateInfo st new mIterations step | step <= mIterations = new <| st
                                   | otherwise = case viewr st of
                                                      (xs :> _) -> new <| xs

