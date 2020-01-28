module Bio.Tools.Sequence.OligoDesigner.Algo
 (designOligsDNA
 ,designOligsAA
 ,getRandomSeed
 ) where

import           Bio.NucleicAcid.Nucleotide.Type            (DNA (..))
import           Bio.Protein.AminoAcid                      (AA)
import qualified Bio.Tools.Sequence.CodonOptimization       as CodonOptimization (optimizeCodonForAA)
import           Bio.Tools.Sequence.OligoDesigner.RnaCofoldOptimizer (minMaxOptimize)
import           Bio.Tools.Sequence.OligoDesigner.Splitter  (split)
import           Bio.Tools.Sequence.OligoDesigner.Types     (OligSet (..), OligSplittingConfig (..),
                                                             OligoDesignerConfig (..))
import           Bio.Tools.Sequence.OligoDesigner.Utils     (buildOligSet)
import           Control.Monad.Except                       (Except, throwError)
import           Control.Monad.IO.Class                     (MonadIO, liftIO)
import           Control.Monad.State                        (State, evalState)
import           System.Random                              (StdGen, getStdGen)

designOligsDNA :: OligSplittingConfig -> [DNA] -> Except String OligSet
designOligsDNA conf sequ =
    case split conf (length sequ) of
        Just splitting -> return $ buildOligSet splitting sequ
        _              -> throwError "Cannot find splitting for parameters"

getRandomSeed :: MonadIO m => m StdGen
getRandomSeed = liftIO getStdGen

designOligsAA :: StdGen -> OligoDesignerConfig -> [AA] -> Except String OligSet
designOligsAA gen conf@(OligoDesignerConfig codonConf _ splittingConf) aa = do
    let dna = CodonOptimization.optimizeCodonForAA codonConf aa
    oligs <- designOligsDNA splittingConf dna
    let state = findOptimal oligs []
    return $ evalState state gen
  where
    findOptimal :: OligSet -> [(OligSet, Double)] -> State StdGen OligSet
    findOptimal oligs _ = minMaxOptimize conf oligs --TODO: do optimization while score isn't stable
