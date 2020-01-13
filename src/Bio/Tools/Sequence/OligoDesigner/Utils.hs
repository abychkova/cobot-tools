module Bio.Tools.Sequence.OligoDesigner.Utils
 (translateDNA
 ,assemble
 ,weightedRandom
 ,randomCodon
 ) where

import Bio.NucleicAcid.Nucleotide.Type (DNA(..))
import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligSet(..))
import Bio.Tools.Sequence.CodonOptimization.Constants (codonFrequencies, ak2Codon)
import Debug.Trace (trace)
import Bio.Protein.AminoAcid.Type (AA)
import Bio.Tools.Sequence.CodonOptimization.Types (Organism)
import Data.Maybe (fromMaybe)
import Data.Map as Map (lookup)
import System.Random (randomRIO)
import Debug.Trace (trace)
import Data.List (sortOn)

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

randomCodon :: Organism -> AA -> IO [DNA]
randomCodon organism aa = weightedRandom codonToWeight where
    codons = fromMaybe [] (Map.lookup aa ak2Codon)
    codonToWeight = map (\codon -> (codon, codonWeight codon)) codons

    codonWeight :: [DNA] -> Double
    codonWeight codon = fromMaybe 0 (Map.lookup codon (codonFrequencies organism))

weightedRandom :: Show a => [(a, Double)] -> IO a
weightedRandom []    = fail "cannot get random for empty array"
weightedRandom items = do
    random <- randomRIO (0, sum $ map snd items)
    return $ getWeightedItem random 0 (sortOn snd items)
  where
    getWeightedItem :: Show a => Double -> Double -> [(a, Double)] -> a
    getWeightedItem randomValue acc ((x, w) : xs) = if acc + w >= randomValue then x else getWeightedItem randomValue (acc + w) xs