module Bio.Tools.Sequence.OligoDesigner.Optimizer.IterativeOptimizer(
    optimize
) where

import Bio.Tools.Sequence.OligoDesigner.Scorer (commonScore)
import Bio.Tools.Sequence.OligoDesigner.Printer (buildStr')
import Bio.Tools.Sequence.OligoDesigner.Optimizer.GCContentOptimizer (gcContentOptimize)
import Bio.Tools.Sequence.OligoDesigner.Optimizer.RNACofoldOptimizer (rnaOptimize)
import System.Random (StdGen)
import Control.Monad.State (State)
import Bio.Tools.Sequence.OligoDesigner.Types (OligsDesignerConfig(..), OligSet)
import Debug.Trace 
import Text.Regex.TDFA (Regex)


optimize :: OligsDesignerConfig -> [Regex] -> OligSet -> State StdGen OligSet
optimize conf regexes oligs = trace (buildStr' conf "" oligs) optimizeIteration 0 [(oligs, commonScore conf oligs)]
  where
    optimizationStep :: OligSet -> State StdGen [(OligSet, Double)]
    optimizationStep oligs =  do
       result <- rnaOptimize conf regexes oligs >>= gcContentOptimize conf regexes
       return $ trace (buildStr' conf "" result) [(result, commonScore conf result)]

    optimizeIteration :: Int -> [(OligSet, Double)] -> State StdGen OligSet
    optimizeIteration iteration results | iteration >= maxOptimizationIteration conf =
                                                trace ("We are reached maximum iteration count. results:" ++ show (map snd results)) $ return $ fst $ last results
                                        | iteration > 2 && isStableScore results =
                                                trace ("We are reached stable score. results:" ++ show (map snd results)) $
                                                return $ fst $ last results
                                        | otherwise =
                                                trace ("iteration № " ++ show iteration) $
                                                optimizationStep (fst $ last results) >>= optimizeIteration (iteration + 1) . (results ++)

isStableScore :: [(OligSet, Double)] -> Bool
isStableScore results = isEqual lastScore prevScore && isEqual lastScore prevPrevScore
  where
    lastScore = snd $ last results
    prevScore = snd $ results !! (length results - 2)
    prevPrevScore = snd $ results !! (length results - 3)
    
isEqual :: Double -> Double -> Bool
isEqual a b = abs(a - b) < 0.00001