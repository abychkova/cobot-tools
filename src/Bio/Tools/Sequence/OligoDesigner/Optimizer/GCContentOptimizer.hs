module Bio.Tools.Sequence.OligoDesigner.Optimizer.GCContentOptimizer where

import Bio.NucleicAcid.Nucleotide (DNA(..))
import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligoDesignerConfig(..), OligSet(..))
import Bio.Tools.Sequence.OligoDesigner.Utils (mutateSlice, assemble, getAANumber, mutate, buildOligSet)
import Bio.Tools.Sequence.OligoDesigner.Scorer (gcContent, gcContentDifference)
import Control.Monad.State (State)
import System.Random (StdGen)
import Bio.Tools.Sequence.CodonOptimization.Types (gcContentDesired, organism)
import Data.List (maximumBy, minimumBy, findIndex)
import Bio.Tools.Sequence.CodonOptimization (CodonOptimizationConfig(..))
import Debug.Trace (trace)

gcContentOptimize :: CodonOptimizationConfig -> OligSet -> State StdGen OligSet
gcContentOptimize (CodonOptimizationConfig organism _ _ _ _ _ _ _ _ _ _ targetGC _)  oligs@(OligSet fwd rvsd splitting) = do
    let allOligs = fwd ++ rvsd
    let maximumOlig = maximumBy gcContentComarator allOligs
    let minimumOlig = minimumBy gcContentComarator allOligs
    let farthestFromTarget =
            if distanceToTarget minimumOlig > distanceToTarget maximumOlig then minimumOlig else maximumOlig

    let indexes = (getAANumber $ start farthestFromTarget, getAANumber $ end farthestFromTarget - 1)
    let dna = assemble oligs

    varSequences <- mutate organism dna indexes
    let variants = map (buildOligSet splitting) varSequences
    return $ minimumBy scoreCmp variants
  where
    
    distanceToTarget :: Olig -> Double
    distanceToTarget (Olig sequ _ _) = abs (gcContent sequ - targetGC)

    scoreCmp :: OligSet -> OligSet -> Ordering
    scoreCmp oligs1 oligs2 = compare score1 score2
      where
        score1 = gcContentDifference oligs1
        score2 = gcContentDifference oligs2

gcContentComarator :: Olig -> Olig -> Ordering
gcContentComarator o1 o2 = compare (gcContent $ sequ o1) (gcContent $ sequ o2)