module Bio.Tools.Sequence.OligoDesigner.Optimizer.GCContentOptimizer where

import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligsDesignerInnerConfig(..), OligSet(..))
import Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils (assemble, getAAIndex, buildOligSet, orderByScore)
import Bio.Tools.Sequence.OligoDesigner.Utils.MutationUtils (mutate)
import Bio.Tools.Sequence.OligoDesigner.Scorer (gcContent, oligsGCContentDifference)
import System.Random (StdGen)
import Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer (filterForbidden)
import Control.Monad.Trans.State.Lazy (StateT)
import Control.Monad.Except (Except)
import Control.Monad.Trans (lift)

gcContentOptimize :: OligsDesignerInnerConfig -> OligSet -> StateT StdGen (Except String) OligSet
gcContentOptimize (OligsDesignerInnerConfig organism targetGC regexes _ _) oligs@(OligSet fwd rvsd splitting) = do
    let allOligs = fwd ++ rvsd
    let (minimumOlig, maximumOlig) = orderByScore allOligs (gcContent . sequDNA)
    let farthestFromTarget = if distanceToTarget minimumOlig > distanceToTarget maximumOlig then minimumOlig else maximumOlig
    startAAIndex <- lift $ getAAIndex $ start farthestFromTarget
    endAAIndex <- lift $ getAAIndex $ end farthestFromTarget - 1
    let dna = assemble oligs
    varSequences <- mutate organism dna (startAAIndex, endAAIndex)
    let filtered =  filterForbidden regexes varSequences
    variants <- lift $ mapM (buildOligSet splitting) filtered
    let (minOligs, _) = orderByScore variants oligsGCContentDifference
    return minOligs
  where
    distanceToTarget :: Olig -> Double
    distanceToTarget (Olig sequ _ _) = abs (gcContent sequ - targetGC)
