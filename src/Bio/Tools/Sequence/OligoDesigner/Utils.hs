module Bio.Tools.Sequence.OligoDesigner.Utils
 (translateDNA
 ,assemble
 ) where

import Bio.NucleicAcid.Nucleotide.Type (DNA(..))
import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligSet(..))
import Debug.Trace (trace)

translateDNA :: [DNA] -> [DNA]
translateDNA = map translate
  where
    translate :: DNA -> DNA
    translate DA = DT
    translate DT = DA
    translate DC = DG
    translate DG = DC

assemble :: OligSet -> [DNA]
assemble (OligSet fwd rvd) = constract fwd rvd 0 [] where
    constract :: [Olig] -> [Olig] -> Int -> [DNA] -> [DNA]
    constract _ [] _ acc = acc
    constract [] _ _ acc = acc
    constract (Olig seq1 startLeft endLeft : xs) (Olig seq2 startRight endRight : ys) prevEnd acc = constract xs ys endRight res
      where
        res = acc ++ drop (prevEnd - startLeft) seq1 ++ translateDNA (drop (endLeft - startRight) seq2)