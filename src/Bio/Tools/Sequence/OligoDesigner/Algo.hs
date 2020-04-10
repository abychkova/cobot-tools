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
                                                                                    OligsDesignerConfig (..), OligsDesignerInnerConfig(..))
import           Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils                         (buildOligSet)
import           Control.Monad.Except                                           (Except, throwError)
import           Control.Monad.IO.Class                                         (MonadIO, liftIO)
import           System.Random                                                  (StdGen, getStdGen)
import           Bio.Tools.Sequence.OligoDesigner.Scorer                        (commonScore)
import Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer(fixForbidden)
import Text.Regex.TDFA (Regex, makeRegex)
import Debug.Trace
import Control.Monad.Trans.State.Lazy (StateT, evalStateT)

designOligsDNA :: OligsSplittingConfig -> [DNA] -> Except String OligSet
designOligsDNA conf dna =
    case split conf (length dna) of
        Just splitting -> return $ buildOligSet splitting dna
        _              -> throwError "Cannot find splitting for parameters"

getRandomSeed :: MonadIO m => m StdGen
getRandomSeed = liftIO getStdGen

designOligsAA :: StdGen -> OligsDesignerConfig -> [AA] -> Except String OligSet
designOligsAA gen conf@(OligsDesignerConfig codonConf _ _ _) aa = do
    let dna = optimizeCodonForAA codonConf aa
    evalStateT (runOptimization conf dna) gen
  where
    runOptimization :: OligsDesignerConfig -> [DNA] -> StateT StdGen (Except String) OligSet
    runOptimization conf@(OligsDesignerConfig _ splittingConf _ _) dna = do
        let innerConf = convertConfig conf
        dnaFixed <- fixForbidden innerConf dna
        case split splittingConf (length dna) of
                Just splitting -> optimize innerConf (buildOligSet splitting dna)
                _              -> throwError "Cannot find splitting for parameters"
                
    convertConfig :: OligsDesignerConfig -> OligsDesignerInnerConfig
    convertConfig conf = OligsDesignerInnerConfig organismType gcTarget regexes optIteration fixForbIteration
      where
        codonConf = codonOptimizationConfig conf
        organismType = organism codonConf
        gcTarget = gcContentDesired codonConf
        regexes = map makeRegex (forbiddenSequence codonConf) :: [Regex]
        optIteration = maxOptimizationIteration conf
        fixForbIteration = maxFixForbiddenIteration conf
        