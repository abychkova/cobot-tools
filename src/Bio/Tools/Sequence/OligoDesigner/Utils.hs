module Bio.Tools.Sequence.OligoDesigner.Utils
 (translate
 ) where

import Bio.NucleicAcid.Nucleotide.Type (DNA(..))

translate :: [DNA] -> [DNA]
translate = map translateDNA
  where
    translateDNA :: DNA -> DNA
    translateDNA DA = DT
    translateDNA DT = DA
    translateDNA DC = DG
    translateDNA DG = DC
