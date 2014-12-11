
module L-BFGS (
              lBFGS 
              ) where

import qualified Data.Vector.Unboxed as U
import Data.Array.Repa as R
import Data.Array.Repa.Algorithms.Matrix (
                        mmultP
                       ,mmultS
                       ,transpose2P
                       ,transpose2S
                       )
import Text.Printf (printf)

-- Internal Modules
import TypesOptimization (
                          Function
                         ,FunGrad
                         ,Gradient
                         ,Matrix
                         ,MaxSteps
                         ,Point
                         ,Tolerance
                          )
import Tools (
             dot
            ,identity
            ,normVec
            ,scalarMatrix
            ,unboxed2Mtx
             )

-- ========================> <===============

lBFGS :: Monad m => Function  -- | Objective funxtion f(X) where X = {x1,x2...xn} 
                 -> FunGrad   -- | Objective function Gradient Df(X) = {df/dx1,..df/dxn} 
                 -> Point     -- | Initial point 
                 -> Matrix    -- | Initial Hessian Matrix
                 -> Tolerance -- | Numerical Tolerance 
                 -> MaxSteps  -- | Maximum allowed steps
                 -> m (Either String Point)
lBFGS f gradF point guessH delta maxSteps = recBFGS point guessH 0
  where recBFGS !xs !hs step =
          if (step > maxSteps)
             then return . Left $ printf "Convergence criterion not met after %d steps\nlast Step:%s\n" step (show xs)
             else do       
              let norma = normVec $ gradF xs        
              if (norma < delta) then return $ Right xs
                                 else do                                      
                                    let d        = 
                                        xsM      = unboxed2Mtx xs
                                        alpha    = wolfeLineSearch f gradF (R.toUnboxed d) xs 2 Nothing
                                        alphaDir = scalarMatrix alpha d 
                                        newXs    = R.toUnboxed .computeUnboxedS $ xsM +^ alphaDir

                                    
                                    
