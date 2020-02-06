module Bio.Tools.Sequence.OligoDesigner.Prettifier(
    prettyDNA
   ,prettySplitting
   ,prettyOlig
   ,prettyOligSet
) where

import Bio.NucleicAcid.Nucleotide (DNA(..))
import Bio.Tools.Sequence.OligoDesigner.Types (OligSplitting(..), Olig(..), OligSet(..))
import Data.List (foldl')


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

mixOligs :: OligSet -> [Olig]
mixOligs (OligSet forward reversed _) = mix forward reversed
  where
    mix :: [a] -> [a] -> [a]
    mix (x:xs) (y:ys) = x : y : mix xs ys
    mix x []          = x
    mix [] y          = y