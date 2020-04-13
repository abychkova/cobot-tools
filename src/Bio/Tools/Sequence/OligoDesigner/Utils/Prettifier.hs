module Bio.Tools.Sequence.OligoDesigner.Utils.Prettifier(
    prettyDNA
   ,prettySplitting
   ,prettyOlig
   ,prettyOligSet
   ,prettyMatrixCell
) where

import Bio.NucleicAcid.Nucleotide (DNA(..))
import Bio.Tools.Sequence.OligoDesigner.Types (OligSplitting(..), Olig(..), OligSet(..), MatrixCell(..))
import Data.List (foldl')
import Data.Matrix (Matrix, toLists)
import Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils (mixOligs)


prettyDNA :: [DNA] -> String
prettyDNA = map prettyOneDNA
  where
    prettyOneDNA :: DNA -> Char
    prettyOneDNA DA = 'A'
    prettyOneDNA DT = 'T'
    prettyOneDNA DC = 'C'
    prettyOneDNA DG = 'G'

prettySplitting :: OligSplitting -> String
prettySplitting (OligSplitting strand5 strand3) = "5': " ++ str5 ++ "\n3': " ++ str3 where
    str5 = fst $ foldl' conc ("", 0) strand5
    str3 = fst $ foldl' conc ("", 0) strand3

    conc :: (String, Int) -> (Int, Int) -> (String, Int)
    conc (res, prev) (x, y) = (res ++ replicate (x - prev) ' ' ++ "(" ++ replicate (y - x) '_' ++ ")", y)

prettyOlig :: Olig -> String
prettyOlig (Olig dna start end) = prettyDNA dna ++ " [" ++ show start ++ ", " ++ show end ++ ")\n"

prettyOligSet :: OligSet -> String
prettyOligSet oligs = concatMap prettyOlig (mixOligs oligs)

prettyMatrixCell :: Matrix MatrixCell -> String
prettyMatrixCell mtx = str
  where
    str = concatMap rowToStr (toLists mtx)

    rowToStr :: [MatrixCell] -> String
    rowToStr cells = concatMap cellToStr cells  ++ "\n"

    cellToStr :: MatrixCell -> String
    cellToStr (MatrixCell _ _ rna) = show rna ++ "\t"