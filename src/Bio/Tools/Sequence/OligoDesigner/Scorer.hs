module Bio.Tools.Sequence.OligoDesigner.Scorer
    ( rnaScore
    , commonScore
    , rnaMatrixScore
    , oligsGCContentDifference
    , gcContent
    , gcContentScoreByOligs
    , gcContentScoreBySequence
    ) where

import           Bio.NucleicAcid.Nucleotide                              (DNA (..))
import           Bio.Tools.Sequence.OligoDesigner.Types                  (MatrixCell (..),
                                                                          Olig (..),
                                                                          OligSet (..),
                                                                          TargetGC)
import           Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils      (assemble)
import           Bio.Tools.Sequence.OligoDesigner.Utils.RNAMatrixBuilder (rnaMatrix)
import           Data.Matrix                                             (Matrix,
                                                                          ncols,
                                                                          nrows,
                                                                          (!))

commonScore :: TargetGC -> OligSet -> Double
commonScore targetGC oligs = scoreValue
  where
    rnaScoreValue = realToFrac $ rnaScore oligs
    oligsGCValue = oligsGCContentDifference oligs
    gcScoreValue = gcContentScoreByOligs oligs targetGC
    scoreValue = rnaScoreValue * gcScoreValue / oligsGCValue

rnaMatrixScore :: Matrix MatrixCell -> Float
rnaMatrixScore oligsMatrix = aboveDiagonalScore - otherScore
  where
    rowsCnt = nrows oligsMatrix
    colsCnt = ncols oligsMatrix
    aboveDiagonalScore = minimum [abs $ rna $ oligsMatrix ! (x , x + 1) | x <- [1 .. rowsCnt - 1]]
    otherScore = maximum [abs $ rna $ oligsMatrix ! (x , y) | x <- [1 .. rowsCnt], y <- [1 .. colsCnt], x <= y && abs (x - y) /= 1]

rnaScore :: OligSet -> Float
rnaScore oligs = rnaMatrixScore $ rnaMatrix oligs

gcContentScoreByOligs :: OligSet -> TargetGC -> Double
gcContentScoreByOligs oligs = gcContentScoreBySequence (assemble oligs)

gcContentScoreBySequence :: [DNA] -> TargetGC -> Double
gcContentScoreBySequence dna targetGC = score where
    deltaGC = abs(targetGC - gcContent dna)
    score = 1 - deltaGC / targetGC -- 1 хорошо и 0 плохо

oligsGCContentDifference :: OligSet -> Double
oligsGCContentDifference (OligSet [] [] _)     = 0
oligsGCContentDifference (OligSet fwrd rvsd _) = maximum gcContents - minimum gcContents
  where
    gcContents = map (gcContent . sequDNA) (fwrd ++ rvsd)

gcContent :: [DNA] -> Double
gcContent sequ =  100 * gc / realToFrac (length sequ)
  where
    gc = realToFrac $ length $ filter (\nk -> nk == DC || nk == DG) sequ
