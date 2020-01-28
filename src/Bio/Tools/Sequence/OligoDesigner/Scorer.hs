module Bio.Tools.Sequence.OligoDesigner.Scorer
 (rnaScore
 ,commonScore
 ,rnaMatrixScore
 ,rnaMatrix
 ,gcContentDifference
 ,gcContent
 ) where

import qualified Bio.Tools.Sequence.CodonOptimization         as CodonOptimization
                                                                                    (score)
import           Bio.Tools.Sequence.OligoDesigner.Types       (MatrixCell (..),
                                                               Olig (..),
                                                               OligSet (..),
                                                               OligoDesignerConfig (..),
                                                               standardTemperature)
import           Bio.Tools.Sequence.OligoDesigner.Utils       (assemble)
import           Bio.Tools.Sequence.ViennaRNA.Internal.Cofold (cofold)
import           Data.Matrix                                  (Matrix (..),
                                                               matrix, (!))
import Bio.NucleicAcid.Nucleotide (DNA(..))
--TODO: test me
commonScore :: OligoDesignerConfig -> OligSet -> Double
commonScore (OligoDesignerConfig codonOptimizationConf balanceFactor _) oligs = scoreValue
  where
    sequ = assemble oligs
    codonScore = CodonOptimization.score codonOptimizationConf sequ
    rnaValue = realToFrac $ rnaMatrixScore $ rnaMatrix oligs
    scoreValue = balanceFactor * codonScore + (1 - balanceFactor) * rnaValue

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
gcContentScore :: [DNA] -> Double -> Double
gcContentScore sequ target = 1 - deltaGC / target 
  where
    deltaGC = abs(target - gcContent sequ)

--TODO: test me
gcContentDifference :: OligSet -> Double
gcContentDifference (OligSet [] [] _)     = 0
gcContentDifference (OligSet fwrd rvsd _) = abs (max - min)
  where
    gcContents = map (gcContent . sequ) (fwrd ++ rvsd)
    max = maximum gcContents
    min = minimum gcContents

--TODO: test me    
gcContent :: [DNA] -> Double
gcContent sequ =  100 * gc / realToFrac (length sequ)
  where
    gc = realToFrac $ length $ filter (\nk -> nk == DC || nk == DG) sequ
