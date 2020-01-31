module Bio.Tools.Sequence.OligoDesigner.Scorer
 (rnaScore
 ,rnaMatrix'
 ,commonScore
 ,rnaMatrixScore
 ,rnaMatrix
 ,gcContentDifference
 ,gcContent
 ,oligsGCContentScore
 ,dnaGCContentScore
 ) where

import qualified Bio.Tools.Sequence.CodonOptimization         as CodonOptimization (score, CodonOptimizationConfig(..))
import           Bio.Tools.Sequence.OligoDesigner.Types       (MatrixCell (..),
                                                               Olig (..),
                                                               OligSet (..),
                                                               OligsDesignerConfig (..),
                                                               standardTemperature)
import           Bio.Tools.Sequence.OligoDesigner.Utils       (assemble)
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
rnaMatrix (OligSet forward reversed _) = matrix rowCount rowCount generator
  where
    allOligs = mix forward reversed
    rowCount = length allOligs

    mix :: [a] -> [a] -> [a]
    mix (x:xs) (y:ys) = x : y : mix xs ys
    mix x []          = x
    mix [] y          = y

    generator :: (Int, Int) -> MatrixCell
    generator (i, j) | i > j     = MatrixCell (Olig "" 0 0 ) (Olig "" 0 0 ) 0 --ignore all above diagonal
                     | otherwise = MatrixCell olig1 olig2 score
      where
        olig1 = allOligs !! (i - 1)
        olig2 = allOligs !! (j - 1)
        score = fst $ cofold standardTemperature (sequ olig1, sequ olig2)

rnaMatrix' :: [[DNA]] -> Matrix Float
rnaMatrix' allOligs = matrix rowCount rowCount generator
  where
    rowCount = length allOligs

    generator :: (Int, Int) -> Float
    generator (i, j) | i > j     = 0 --ignore all above diagonal
                     | otherwise = score
      where
        olig1 = allOligs !! (i - 1)
        olig2 = allOligs !! (j - 1)
        score = fst $ cofold standardTemperature (olig1, olig2)

rnaMatrixScore :: Matrix MatrixCell -> Float
rnaMatrixScore oligsMatrix = aboveDiagonalScore - otherScore
  where
    rowsCnt = nrows oligsMatrix
    colsCnt = ncols oligsMatrix
    aboveDiagonalScore = minimum [rna $ oligsMatrix ! (x , y) | x <- [1 .. rowsCnt], y <- [1 .. colsCnt], abs (x - y) == 1]
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
