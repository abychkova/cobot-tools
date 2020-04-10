module Bio.Tools.Sequence.OligoDesigner.Optimizer.IterativeOptimizer(
    optimize
) where

import Bio.Tools.Sequence.OligoDesigner.Scorer (commonScore)
import Bio.Tools.Sequence.OligoDesigner.Utils.Printer (buildStr')
import Bio.Tools.Sequence.OligoDesigner.Optimizer.GCContentOptimizer (gcContentOptimize)
import Bio.Tools.Sequence.OligoDesigner.Optimizer.RNACofoldOptimizer (rnaOptimize)
import System.Random (StdGen)
import Control.Monad.State (State)
import Bio.Tools.Sequence.OligoDesigner.Types (OligsDesignerInnerConfig(..), OligSet)
import Debug.Trace 
import Text.Regex.TDFA (Regex)
import Control.Monad.Trans.State.Lazy (StateT)
import Control.Monad.Except (Except)


optimize :: OligsDesignerInnerConfig -> OligSet -> StateT StdGen (Except String) OligSet
optimize conf@(OligsDesignerInnerConfig _ targetGC _ maxIteration _) oligs = 
    optimizeIteration 0 [(oligs, commonScore targetGC oligs)]
  where
    optimizationStep :: OligSet -> StateT StdGen (Except String) [(OligSet, Double)]
    optimizationStep oligs =  do
       result <- rnaOptimize conf oligs >>= gcContentOptimize conf
       return [(result, commonScore targetGC result)]

    optimizeIteration :: Int -> [(OligSet, Double)] -> StateT StdGen (Except String) OligSet
    optimizeIteration iteration results | iteration >= maxIteration =
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