module Bio.Tools.Sequence.OligoDesigner.Utils.RNAMatrixBuilder(
    rnaMatrix,
    rebuildMatrix
) where


import Data.Matrix (Matrix, matrix, Matrix(..), (!))
import Bio.Tools.Sequence.OligoDesigner.Types (OligSet, Olig(..), MatrixCell(..), OligLight(..), emptyMatrixCell, standardTemperature, forward, reversed)
import Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils (mixOligs)
import Bio.Tools.Sequence.OligoDesigner.Utils.Prettifier (prettyDNA)
import Bio.Tools.Sequence.ViennaRNA.Cofold (cofold)
import Debug.Trace (trace)

rnaMatrix :: OligSet -> Matrix MatrixCell
rnaMatrix oligs = matrix rowCount rowCount generator
  where
    allOligs = mixOligs oligs
    rowCount = length allOligs

    generator :: (Int, Int) -> MatrixCell
    generator (i, j) | i > j     = emptyMatrixCell --ignore all above diagonal
                     | otherwise = MatrixCell (OligLight (prettyDNA dna1) olig1) (OligLight (prettyDNA dna2) olig2) rna
      where
        olig1@(Olig dna1 _ _) = allOligs !! (i - 1)
        olig2@(Olig dna2 _ _) = allOligs !! (j - 1)
        rna = fst $ cofold standardTemperature (dna1, dna2)

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
        olig1@(Olig dna1 start1 end1) = allOligs !! (i - 1)
        olig2@(Olig dna2 start2 end2) = allOligs !! (j - 1)

        olig1L = OligLight (prettyDNA dna1) olig1
        olig2L = OligLight (prettyDNA dna2) olig2

        newMatrixCell = if oldDna1 == dna1 && oldDna2 == dna2
            then oldCell
            else MatrixCell olig1L olig2L (fst $ cofold standardTemperature (dna1, dna2))
