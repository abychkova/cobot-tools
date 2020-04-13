module Bio.Tools.Sequence.OligoDesigner.Optimizer.GCContentOptimizer where

import           Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer      (filterForbidden)
import           Bio.Tools.Sequence.OligoDesigner.Scorer              (gcContent,
                                                                       oligsGCContentDifference)
import           Bio.Tools.Sequence.OligoDesigner.Types               (Olig (..),
                                                                       OligSet (..),
                                                                       OligsDesignerInnerConfig (..)
                                                                       )
import           Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils   (assemble, buildOligSet,
                                                                       getAAIndex,
                                                                       orderByScore)
import           Bio.Tools.Sequence.OligoDesigner.Utils.MutationUtils (mutate)
import           Control.Monad.Except                                 (Except)
import           Control.Monad.Trans                                  (lift)
import           Control.Monad.Trans.State.Lazy                       (StateT)
import           System.Random                                        (StdGen)

-- | 'gcContentOptimize' function does optimization for oligs gc-contents.
-- The main idea is to reduce the difference between oligs gc-content and make it closer to the target
gcContentOptimize :: OligsDesignerInnerConfig    -- ^ configuration data (used 'organism', 'targetGC' and 'regexes')
        -> OligSet                               -- ^ oligs set
        -> StateT StdGen (Except String) OligSet -- ^ result of optimization is optimized oligs set with the best score or error string
gcContentOptimize (OligsDesignerInnerConfig organism targetGC regexes _ _) oligs@(OligSet fwd rvsd splitting) = do
    let allOligs = fwd ++ rvsd
    let (minimumOlig, maximumOlig) = orderByScore allOligs (gcContent . sequDNA)
    let farthestFromTarget =
            if distanceToTarget minimumOlig > distanceToTarget maximumOlig
                then minimumOlig
                else maximumOlig
    startAAIndex <- lift $ getAAIndex $ start farthestFromTarget
    endAAIndex <- lift $ getAAIndex $ end farthestFromTarget - 1
    let dna = assemble oligs
    varSequences <- mutate organism dna (startAAIndex, endAAIndex)
    let filtered = filterForbidden regexes varSequences
    variants <- lift $ mapM (buildOligSet splitting) filtered
    let (minOligs, _) = orderByScore variants oligsGCContentDifference
    return minOligs
  where
    distanceToTarget :: Olig -> Double
    distanceToTarget (Olig sequ _ _) = abs (gcContent sequ - targetGC)