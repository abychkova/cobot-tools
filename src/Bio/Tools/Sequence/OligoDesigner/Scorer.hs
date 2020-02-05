module Bio.Tools.Sequence.OligoDesigner.Scorer
 (rnaScore
 ,commonScore
 ,rnaMatrixScore
 ,rnaMatrix
 ,gcContentDifference
 ,gcContent
 ,oligsGCContentScore
 ,dnaGCContentScore
 ,rebuildMatrix
 ) where

import qualified Bio.Tools.Sequence.CodonOptimization         as CodonOptimization (score, CodonOptimizationConfig(..))
import           Bio.Tools.Sequence.OligoDesigner.Types       (MatrixCell (..),
                                                               Olig (..),
                                                               OligSet (..),
                                                               OligsDesignerConfig (..),
                                                               standardTemperature)
import           Bio.Tools.Sequence.OligoDesigner.Utils       (assemble, mixOligs)
import           Bio.Tools.Sequence.ViennaRNA.Internal.Cofold (cofold)
import           Data.Matrix                                  (Matrix (..),
                                                               matrix, (!))
import Bio.NucleicAcid.Nucleotide (DNA(..))
import Debug.Trace (trace)
import Bio.Tools.Sequence.CodonOptimization.Types (gcContentDesired)

commonScore :: OligsDesignerConfig -> OligSet -> Double
commonScore (OligsDesignerConfig codonConf _ rnaF oligsGCF gcF _) oligs = scoreValue
  where
    rnaScoreValue = realToFrac $ rnaScore oligs
    oligsGCValue = gcContentDifference oligs
    gcScoreValue = oligsGCContentScore oligs (gcContentDesired codonConf)
    scoreValue = rnaF * rnaScoreValue + oligsGCF * oligsGCValue + gcF * gcScoreValue

rnaMatrix :: OligSet -> Matrix MatrixCell
rnaMatrix oligs = matrix rowCount rowCount generator
  where
    allOligs = mixOligs oligs
    rowCount = length allOligs

    generator :: (Int, Int) -> MatrixCell
    generator (i, j) | i > j     = MatrixCell (Olig "" 0 0 ) (Olig "" 0 0 ) 0 --ignore all above diagonal
                     | otherwise = MatrixCell olig1 olig2 score
      where
        olig1 = allOligs !! (i - 1)
        olig2 = allOligs !! (j - 1)
        score = fst $ cofold standardTemperature (sequ olig1, sequ olig2)

rebuildMatrix :: Matrix MatrixCell -> OligSet -> Matrix MatrixCell
rebuildMatrix oldMatrix oligs = newMatrix
  where
    newMatrix = matrix rowCount rowCount generator
    allOligs = mixOligs oligs
    rowCount = length allOligs

    generator :: (Int, Int) -> MatrixCell
    generator (i, j) | i > j     = MatrixCell (Olig "" 0 0 ) (Olig "" 0 0 ) 0
                     | otherwise = newMatrixCell
      where
        oldCell@(MatrixCell oldOlig1 oldOlig2 oldRna) = oldMatrix ! (i, j)
        olig1 = allOligs !! (i - 1)
        olig2 = allOligs !! (j - 1)
        newMatrixCell = if oldOlig1 == olig1 && oldOlig2 == olig2
            then oldCell
            else MatrixCell olig1 olig2 (fst $ cofold standardTemperature (sequ olig1, sequ olig2))

rnaMatrixScore :: Matrix MatrixCell -> Float
rnaMatrixScore oligsMatrix = aboveDiagonalScore - otherScore
  where
    rowsCnt = nrows oligsMatrix
    colsCnt = ncols oligsMatrix
    aboveDiagonalScore = minimum [rna $ oligsMatrix ! (x , x + 1) | x <- [1 .. rowsCnt - 1]]
    otherScore         = maximum [rna $ oligsMatrix ! (x , y) | x <- [1 .. rowsCnt], y <- [1 .. colsCnt], x <= y && abs (x - y) /= 1]

rnaScore :: OligSet -> Float
rnaScore oligs = rnaMatrixScore $ rnaMatrix oligs

--TODO: test me
oligsGCContentScore :: OligSet -> Double -> Double
oligsGCContentScore oligs = dnaGCContentScore (assemble oligs)

dnaGCContentScore :: [DNA] -> Double -> Double
dnaGCContentScore dna target = score where
    deltaGC = abs(target - gcContent dna)
    score = 1 - deltaGC / target

gcContentDifference :: OligSet -> Double
gcContentDifference (OligSet [] [] _)     = 0
gcContentDifference (OligSet fwrd rvsd _) = maximum gcContents - minimum gcContents
  where
    gcContents = map (gcContent . sequ) (fwrd ++ rvsd)

gcContent :: [DNA] -> Double
gcContent sequ =  100 * gc / realToFrac (length sequ)
  where
    gc = realToFrac $ length $ filter (\nk -> nk == DC || nk == DG) sequ
