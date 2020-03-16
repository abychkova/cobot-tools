module Bio.Tools.Sequence.OligoDesigner.RNAMatrixBuilder(
    rnaMatrix,
    rebuildMatrix
) where


import Data.Matrix (Matrix, matrix, Matrix(..), (!))
import Bio.Tools.Sequence.OligoDesigner.Types (OligSet, Olig(..), MatrixCell(..), OligLight(..), emptyMatrixCell, standardTemperature)
import Bio.Tools.Sequence.OligoDesigner.Utils (mixOligs)
import Bio.Tools.Sequence.OligoDesigner.Prettifier (prettyDNA)
import Bio.Tools.Sequence.ViennaRNA.Cofold (cofold)

rnaMatrix :: OligSet -> Matrix MatrixCell
rnaMatrix oligs = matrix rowCount rowCount generator
  where
    allOligs = mixOligs oligs
    rowCount = length allOligs

    generator :: (Int, Int) -> MatrixCell
    generator (i, j) | i > j     = emptyMatrixCell --ignore all above diagonal
                     | otherwise = MatrixCell (OligLight (prettyDNA dna1) olig1) (OligLight (prettyDNA dna2) olig2) score
      where
        olig1@(Olig dna1 _ _) = allOligs !! (i - 1)
        olig2@(Olig dna2 _ _) = allOligs !! (j - 1)
        score = fst $ cofold standardTemperature (dna1, dna2)

rebuildMatrix :: Matrix MatrixCell -> OligSet -> Matrix MatrixCell
rebuildMatrix oldMatrix oligs = newMatrix
  where
    newMatrix = matrix rowCount rowCount generator
    allOligs = mixOligs oligs
    rowCount = length allOligs

    generator :: (Int, Int) -> MatrixCell
    generator (i, j) | i > j     = emptyMatrixCell
                     | otherwise = newMatrixCell
      where
        oldCell@(MatrixCell (OligLight _ (Olig oldDna1 _ _)) (OligLight _ (Olig oldDna2 _ _)) _) = oldMatrix ! (i, j)
--        oldCell@(MatrixCell (Olig oldDna1 _ _) (Olig oldDna2 _ _) _) = oldMatrix ! (i, j)
        olig1@(Olig dna1 start1 end1) = allOligs !! (i - 1)
        olig2@(Olig dna2 start2 end2) = allOligs !! (j - 1)

        olig1L = OligLight (prettyDNA dna1) olig1
        olig2L = OligLight (prettyDNA dna2) olig2

        --FIXME: слишком страшно с переводом в строку. это реально нужно для скорости?
        newMatrixCell = if oldDna1 == dna1 && oldDna2 == dna2
            then oldCell
--            else MatrixCell olig1 olig2 (fst $ cofold standardTemperature (dna1, dna2))
            else MatrixCell olig1L olig2L (fst $ cofold standardTemperature (dna1, dna2))
