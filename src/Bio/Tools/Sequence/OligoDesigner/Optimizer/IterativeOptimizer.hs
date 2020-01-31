module Bio.Tools.Sequence.OligoDesigner.Optimizer.IterativeOptimizer(
    optimize
) where

import Bio.Tools.Sequence.OligoDesigner.Scorer (commonScore)
import Bio.Tools.Sequence.OligoDesigner.Optimizer.GCContentOptimizer (gcContentOptimize)
import Bio.Tools.Sequence.OligoDesigner.Optimizer.RNACofoldOptimizer (rnaOptimize)
import System.Random (StdGen)
import Control.Monad.State (State)
import Bio.Tools.Sequence.OligoDesigner.Types (OligsDesignerConfig(..), OligSet)
import Debug.Trace (trace)


optimize :: OligsDesignerConfig -> OligSet -> State StdGen OligSet
optimize conf oligs = optimizeIteration 0 [(oligs, commonScore conf oligs)]
  where
    optimizationStep :: OligSet -> State StdGen [(OligSet, Double)]
    optimizationStep oligs =  do
       result1 <- rnaOptimize conf oligs
       result2 <- gcContentOptimize conf result1
       return [(result1, commonScore conf result1), (result2, commonScore conf result2)]

    optimizeIteration :: Int -> [(OligSet, Double)] -> State StdGen OligSet
    optimizeIteration iteration results | iteration >= maxOptimizationIteration conf =
                                                trace "We are recached maximum iteration count" $ return $ fst $ last results
                                        | iteration > 2 && isStableScore results    =
                                                trace "We are recached stable score" $ return $ fst $ last results
                                        | otherwise =
                                                trace ("iteration № " ++ show iteration) $
                                                optimizationStep (fst $ last results) >>= optimizeIteration (iteration + 1) . (results ++)

isStableScore :: [(OligSet, Double)] -> Bool
isStableScore results = lastScore == prevScore && lastScore == prevPrevScore
  where
    lastScore = snd $ last results
    prevScore = snd $ results !! (length results - 3)
    prevPrevScore = snd $ results !! (length results - 5)