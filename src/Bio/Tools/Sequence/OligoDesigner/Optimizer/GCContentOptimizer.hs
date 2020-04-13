module Bio.Tools.Sequence.OligoDesigner.Optimizer.GCContentOptimizer where

import Bio.NucleicAcid.Nucleotide (DNA(..))
import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligsDesignerInnerConfig(..), OligSet(..), TargetGC)
import Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils (assemble, getAAIndex, buildOligSet, orderByScore)
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
import Control.Monad.Trans.State.Lazy (StateT)
import Control.Monad.Except (Except)
import Control.Monad.Trans (lift)

gcContentOptimize :: OligsDesignerInnerConfig -> OligSet -> StateT StdGen (Except String) OligSet
gcContentOptimize (OligsDesignerInnerConfig organism targetGC regexes _ _) oligs@(OligSet fwd rvsd splitting) = do
    let allOligs = fwd ++ rvsd
    let (minimumOlig, maximumOlig) = orderByScore allOligs (gcContent . sequDNA)
    let farthestFromTarget =
            if distanceToTarget targetGC minimumOlig > distanceToTarget targetGC maximumOlig then minimumOlig else maximumOlig

    startAAIndex <- lift $ getAAIndex $ start farthestFromTarget
    endAAIndex <- lift $ getAAIndex $ end farthestFromTarget - 1
    let dna = assemble oligs
    varSequences <- mutate organism dna (startAAIndex, endAAIndex)
    let filtered =  filterForbidden regexes varSequences
    variants <- lift $ mapM (buildOligSet splitting) filtered
    let (min, _) = orderByScore variants oligsGCContentDifference
    return min
  where
    distanceToTarget :: TargetGC -> Olig -> Double
    distanceToTarget targetGC (Olig sequ _ _) = abs (gcContent sequ - targetGC)
