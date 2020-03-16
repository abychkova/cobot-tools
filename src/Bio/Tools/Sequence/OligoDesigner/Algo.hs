module Bio.Tools.Sequence.OligoDesigner.Algo
 (designOligsDNA
 ,designOligsAA
 ,getRandomSeed
 ) where

import           Bio.NucleicAcid.Nucleotide.Type                                (DNA (..))
import           Bio.Protein.AminoAcid                                          (AA)
import           Bio.Tools.Sequence.CodonOptimization                           as CodonOptimization (optimizeCodonForAA, CodonOptimizationConfig(..))
import           Bio.Tools.Sequence.OligoDesigner.Optimizer.IterativeOptimizer (optimize)
import           Bio.Tools.Sequence.OligoDesigner.Optimizer.GCContentOptimizer (gcContentOptimize)
import           Bio.Tools.Sequence.OligoDesigner.Splitter                      (split)
import           Bio.Tools.Sequence.OligoDesigner.Types                         (OligSet (..), OligsSplittingConfig (..),
                                                                                    OligsDesignerConfig (..))
import           Bio.Tools.Sequence.OligoDesigner.Utils                         (buildOligSet)
import           Control.Monad.Except                                           (Except, throwError)
import           Control.Monad.IO.Class                                         (MonadIO, liftIO)
import           Control.Monad.State                                            (State, evalState)
import           System.Random                                                  (StdGen, getStdGen)
import           Bio.Tools.Sequence.OligoDesigner.Scorer                        (commonScore)
import Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer(fixForbidden)
import Text.Regex.TDFA (Regex, makeRegex)
import Debug.Trace

designOligsDNA :: OligsDesignerConfig -> [DNA] -> Except String OligSet
designOligsDNA (OligsDesignerConfig _ conf _ _ _ _ _) dna =
    case split conf (length dna) of
        Just splitting -> return $ buildOligSet splitting dna
        _              -> throwError "Cannot find splitting for parameters"

getRandomSeed :: MonadIO m => m StdGen
getRandomSeed = liftIO getStdGen

designOligsAA :: StdGen -> OligsDesignerConfig -> [AA] -> Except String OligSet
designOligsAA gen conf@(OligsDesignerConfig codonConf _ _ _ _ _ _) aa = do
    let dna = optimizeCodonForAA codonConf aa
    let regexes = map makeRegex (forbiddenSequence codonConf) :: [Regex]
    dnaFixed <- fixForbidden gen conf regexes dna
    oligs    <- designOligsDNA conf dnaFixed
    return $ evalState (optimize conf regexes oligs) gen