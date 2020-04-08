module Bio.Tools.Sequence.OligoDesigner.Scorer
 (rnaScore
 ,commonScore
 ,rnaMatrixScore
 ,oligsGCContentDifference
 ,gcContent
 ,gcContentScoreByOligs
 ,gcContentScoreBySequence
 ) where

import qualified Bio.Tools.Sequence.CodonOptimization         as CodonOptimization (score, CodonOptimizationConfig(..))
import           Bio.Tools.Sequence.OligoDesigner.Types       (MatrixCell (..),
                                                               Olig (..),
                                                               OligLight (..),
                                                               OligSet (..),
                                                               OligsDesignerConfig (..),
                                                               standardTemperature, emptyMatrixCell)
import           Bio.Tools.Sequence.OligoDesigner.Utils       (assemble, mixOligs)
import           Bio.Tools.Sequence.OligoDesigner.RNAMatrixBuilder       (rnaMatrix)
import           Bio.Tools.Sequence.ViennaRNA.Internal.Cofold (cofold)
import Bio.Tools.Sequence.OligoDesigner.Prettifier (prettyDNA, prettyMatrixCell)
import Bio.NucleicAcid.Nucleotide (DNA(..))
import Debug.Trace
import Bio.Tools.Sequence.CodonOptimization.Types (gcContentDesired)
import Data.Text (pack)
import Data.Matrix (Matrix, nrows, ncols, (!), prettyMatrix)

commonScore :: OligsDesignerConfig -> OligSet -> Double
commonScore (OligsDesignerConfig codonConf _ _ _) oligs = scoreValue
  where
    rnaScoreValue = realToFrac $ rnaScore oligs
    oligsGCValue = oligsGCContentDifference oligs
    gcScoreValue = gcContentScoreByOligs oligs (gcContentDesired codonConf)
    scoreValue = rnaScoreValue * gcScoreValue / oligsGCValue
    
rnaMatrixScore :: Matrix MatrixCell -> Float
rnaMatrixScore oligsMatrix = aboveDiagonalScore - otherScore
  where
    rowsCnt = nrows oligsMatrix
    colsCnt = ncols oligsMatrix
    aboveDiagonalScore = minimum [abs $ rna $ oligsMatrix ! (x , x + 1) | x <- [1 .. rowsCnt - 1]]
    otherScore         = maximum [abs $ rna $ oligsMatrix ! (x , y) | x <- [1 .. rowsCnt], y <- [1 .. colsCnt], x <= y && abs (x - y) /= 1]

rnaScore :: OligSet -> Float
rnaScore oligs = rnaMatrixScore $ rnaMatrix oligs

--TODO: test me
gcContentScoreByOligs :: OligSet -> Double -> Double
gcContentScoreByOligs oligs = gcContentScoreBySequence (assemble oligs)

gcContentScoreBySequence :: [DNA] -> Double -> Double
gcContentScoreBySequence dna target = score where
    deltaGC = abs(target - gcContent dna)
    score = 1 - deltaGC / target -- 1 хорошо и 0 плохо

oligsGCContentDifference :: OligSet -> Double
oligsGCContentDifference (OligSet [] [] _)     = 0
oligsGCContentDifference (OligSet fwrd rvsd _) = maximum gcContents - minimum gcContents
  where
    gcContents = map (gcContent . sequDNA) (fwrd ++ rvsd)

gcContent :: [DNA] -> Double
gcContent sequ =  100 * gc / realToFrac (length sequ)
  where
    gc = realToFrac $ length $ filter (\nk -> nk == DC || nk == DG) sequ
