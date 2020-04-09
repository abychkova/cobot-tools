module Bio.Tools.Sequence.OligoDesigner.Optimizer.GCContentOptimizer where

import Bio.NucleicAcid.Nucleotide (DNA(..))
import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligsDesignerInnerConfig(..), OligSet(..), TargetGC)
import Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils (assemble, getAAIndex, buildOligSet, compareBySecond)
import Bio.Tools.Sequence.OligoDesigner.Utils.MutationUtils (mutateSlice, mutate)
import Bio.Tools.Sequence.OligoDesigner.Scorer (gcContent, oligsGCContentDifference)
import Control.Monad.State (State)
import System.Random (StdGen)
import Bio.Tools.Sequence.CodonOptimization.Types (gcContentDesired, organism)
import Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer (filterForbidden)
import Data.List (maximumBy, minimumBy, findIndex)
import Bio.Tools.Sequence.CodonOptimization (CodonOptimizationConfig(..))
import Debug.Trace
import Text.Regex.TDFA (Regex)
import Bio.Tools.Sequence.OligoDesigner.Utils.Prettifier (prettyOlig, prettyDNA)

gcContentOptimize :: OligsDesignerInnerConfig -> OligSet -> State StdGen OligSet
gcContentOptimize (OligsDesignerInnerConfig organism targetGC regexes _ _) oligs@(OligSet fwd rvsd splitting) = do
    let allOligs = fwd ++ rvsd
    let oligsToGc = map (\o -> (o, gcContent $ sequDNA o)) allOligs
    let maximumOlig = fst $ maximumBy compareBySecond oligsToGc
    let minimumOlig = fst $ minimumBy compareBySecond oligsToGc
    let farthestFromTarget =
            if distanceToTarget targetGC minimumOlig > distanceToTarget targetGC maximumOlig then minimumOlig else maximumOlig

    let indexesToMutate = (getAAIndex $ start farthestFromTarget, getAAIndex $ end farthestFromTarget - 1)
    let dna = assemble oligs
    varSequences <- mutate organism dna indexesToMutate
    let filtered =  filterForbidden regexes varSequences
    let variants = map (buildOligSet splitting) filtered
    let varsToScore = map (\vars -> (vars, oligsGCContentDifference vars)) variants
    return $ fst $ minimumBy compareBySecond varsToScore
  where
    distanceToTarget :: TargetGC -> Olig -> Double
    distanceToTarget targetGC (Olig sequ _ _) = abs (gcContent sequ - targetGC)