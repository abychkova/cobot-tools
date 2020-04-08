module Bio.Tools.Sequence.OligoDesigner.Optimizer.GCContentOptimizer where

import Bio.NucleicAcid.Nucleotide (DNA(..))
import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligsDesignerConfig(..), OligSet(..))
import Bio.Tools.Sequence.OligoDesigner.Utils (mutateSlice, assemble, getAAIndex, mutate, buildOligSet, compareBySecond)
import Bio.Tools.Sequence.OligoDesigner.Scorer (gcContent, oligsGCContentDifference)
import Control.Monad.State (State)
import System.Random (StdGen)
import Bio.Tools.Sequence.CodonOptimization.Types (gcContentDesired, organism)
import Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer (filterForbidden)
import Data.List (maximumBy, minimumBy, findIndex)
import Bio.Tools.Sequence.CodonOptimization (CodonOptimizationConfig(..))
import Debug.Trace
import Text.Regex.TDFA (Regex)
import Bio.Tools.Sequence.OligoDesigner.Prettifier (prettyOlig, prettyDNA)

gcContentOptimize :: OligsDesignerConfig -> [Regex] -> OligSet -> State StdGen OligSet
gcContentOptimize conf regexes oligs@(OligSet fwd rvsd splitting) = do
    let targetGC = gcContentDesired $ codonOptimizationConfig conf
    let organismType = organism $ codonOptimizationConfig conf
    let allOligs = fwd ++ rvsd
    let oligsToGc = map (\o -> (o, gcContent $ sequDNA o)) allOligs
    let maximumOlig = fst $ maximumBy compareBySecond oligsToGc
    let minimumOlig = fst $ minimumBy compareBySecond oligsToGc
    let farthestFromTarget =
            if distanceToTarget targetGC minimumOlig > distanceToTarget targetGC maximumOlig then minimumOlig else maximumOlig

    let indexesToMutate = (getAAIndex $ start farthestFromTarget, getAAIndex $ end farthestFromTarget - 1)
    let dna = assemble oligs
    varSequences <- mutate organismType dna indexesToMutate
    let filtered =  filterForbidden regexes varSequences
    let variants = map (buildOligSet splitting) filtered
    let varsToScore = map (\vars -> (vars, oligsGCContentDifference vars)) variants
    return $ fst $ minimumBy compareBySecond varsToScore
  where
    distanceToTarget :: Double -> Olig -> Double
    distanceToTarget targetGC (Olig sequ _ _) = abs (gcContent sequ - targetGC)