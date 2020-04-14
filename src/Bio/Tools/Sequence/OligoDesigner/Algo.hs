module Bio.Tools.Sequence.OligoDesigner.Algo
    ( designOligsDNA
    , designOligsAA
    , getRandomSeed
    ) where

import Control.Monad.Except           (Except, throwError)
import Control.Monad.IO.Class         (MonadIO, liftIO)
import Control.Monad.Trans            (lift)
import Control.Monad.Trans.State.Lazy (StateT, evalStateT)
import System.Random                  (StdGen, getStdGen)
import Text.Regex.TDFA                (Regex, makeRegex)
                                                                                
import Bio.NucleicAcid.Nucleotide.Type (DNA (..))
import Bio.Protein.AminoAcid           (AA)

import Bio.Tools.Sequence.CodonOptimization                          (CodonOptimizationConfig (..),
                                                                      optimizeCodonForAA)
import Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer               (fixForbidden)
import Bio.Tools.Sequence.OligoDesigner.Optimizer.IterativeOptimizer (optimize)
import Bio.Tools.Sequence.OligoDesigner.Splitter                     (split)
import Bio.Tools.Sequence.OligoDesigner.Types                        (OligSet (..), OligSplitting,
                                                                      OligsDesignerConfig (..),
                                                                      OligsDesignerInnerConfig (..),
                                                                      OligsSplittingConfig (..))
import Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils            (buildOligSet)


-- | function does splitting DNA-sequence to oligs according to splitting config.
designOligsDNA :: OligsSplittingConfig -- ^ Splitting configuration
               -> [DNA]                   -- ^ DNA-sequence
               -> Except String OligSet   -- ^ Result is set of oligs or error string
designOligsDNA conf dna =
    case split conf (length dna) of
        Just splitting -> buildOligSet splitting dna
        _              -> throwError "Cannot find splitting for parameters"

-- | function gets random generator. You can use it as parameter for 'designOligsAA' function.
getRandomSeed :: MonadIO m => m StdGen
getRandomSeed = liftIO getStdGen

-- | function does splitting and optimization AA-sequence to oligs according config.
designOligsAA :: StdGen          -- ^ Random generator
              -> OligsDesignerConfig   -- ^ configuration for splitting and for codon-optimization
              -> [AA]                  -- ^ AA-sequence
              -> Except String OligSet -- ^ Result is set of oligs with the best score or error string
designOligsAA gen conf@(OligsDesignerConfig codonConf splittingConf _ _) aa = do
    let dna = optimizeCodonForAA codonConf aa
    case split splittingConf (length dna) of
            Just splitting -> evalStateT (runOptimization dna splitting) gen
            _              -> throwError "Cannot find splitting for parameters"

  where
    runOptimization :: [DNA] -> OligSplitting -> StateT StdGen (Except String) OligSet
    runOptimization dna splitting = do
        let innerConf = convertConfig
        dnaFixed <- fixForbidden innerConf dna
        oligs <- lift $ buildOligSet splitting dnaFixed
        optimize innerConf oligs

    convertConfig :: OligsDesignerInnerConfig
    convertConfig = OligsDesignerInnerConfig organismType gcTarget regexes optIteration fixForbIteration
      where
        organismType = organism codonConf
        gcTarget = gcContentDesired codonConf
        regexes = map makeRegex (forbiddenSequence codonConf) :: [Regex]
        optIteration = maxOptimizationIteration conf
        fixForbIteration = maxFixForbiddenIteration conf

