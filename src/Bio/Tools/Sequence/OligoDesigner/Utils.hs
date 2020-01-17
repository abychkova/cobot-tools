module Bio.Tools.Sequence.OligoDesigner.Utils
 (assemble
 ,weightedRandom
 ,randomCodon
 ,buildOligSet
 ,oneMutation
 ,slice
 ) where

import           Bio.NucleicAcid.Nucleotide.Type                (DNA (..), cNA)
import           Bio.Protein.AminoAcid.Type                     (AA)
import           Bio.Tools.Sequence.CodonOptimization.Constants (ak2Codon,
                                                                 codon2ak,
                                                                 codonFrequencies)
import           Bio.Tools.Sequence.CodonOptimization.Types     (Organism)
import           Bio.Tools.Sequence.OligoDesigner.Types         (Codon,
                                                                 Olig (..),
                                                                 OligBounds,
                                                                 OligSet (..),
                                                                 OligSplitting (..))
import           Control.Monad.State                            (State, get,
                                                                 put)
import           Data.List                                      (sortOn)
import           Data.Map                                       as Map (lookup)
import           Data.Maybe                                     (fromMaybe)
import           System.Random                                  (StdGen,
                                                                 randomR)

assemble :: OligSet -> [DNA]
assemble (OligSet fwd rvd _) = constract fwd rvd 0 [] where
    constract :: [Olig] -> [Olig] -> Int -> [DNA] -> [DNA]
    constract _ [] _ acc = acc
    constract [] _ _ acc = acc
    constract (Olig seq1 startLeft endLeft : xs) (Olig seq2 startRight endRight : ys) prevEnd acc = constract xs ys endRight res
      where
        res = acc ++ drop (prevEnd - startLeft) seq1 ++ map cNA (drop (endLeft - startRight) seq2)

--TODO: test me
--TODO: correct exception instead
oneMutation :: Organism -> Codon -> State StdGen [DNA]
oneMutation organism codon = do
    let aa = fromMaybe (error ("cannot find aa for codon " ++ show codon)) (Map.lookup codon codon2ak)
    newCodon <- randomCodon organism aa
    if newCodon == codon
        then oneMutation organism codon
        else return newCodon

--TODO: test me
randomCodon :: Organism -> AA -> State StdGen Codon
randomCodon organism aa = do
    let codons = fromMaybe [] (Map.lookup aa ak2Codon)
    let codonToWeight = map (\codon -> (codon, codonWeight codon)) codons
    weightedRandom codonToWeight
  where
    codonWeight :: [DNA] -> Double
    codonWeight codon = fromMaybe 0 (Map.lookup codon (codonFrequencies organism))

weightedRandom :: [(a, Double)] -> State StdGen a
weightedRandom []    = error "cannot get random for empty array" --FIXME: normal exception instead
weightedRandom items = do
    gen <- get
    let (random, newGen) = randomR (0, sum $ map snd items) gen
    put newGen
    let value = getWeightedItem random 0 (sortOn snd items)
    return value
  where
    getWeightedItem :: Double -> Double -> [(a, Double)] -> a
    getWeightedItem _ _ []       = error "some strange situation"
    getWeightedItem _ _ [(x, _)] = x
    getWeightedItem randomValue acc ((x, w) : xs) | acc + w >= randomValue = x
                                                  | otherwise = getWeightedItem randomValue (acc + w) xs

buildOligSet :: OligSplitting -> [DNA] -> OligSet
buildOligSet splitting sequ = OligSet strand5' strand3' splitting
  where
    strand5' = map (buildOlig sequ) (strand5 splitting)
    strand3' = map (buildOlig $ map cNA sequ) (strand3 splitting)

    buildOlig :: [DNA] -> OligBounds -> Olig
    buildOlig dna (start, end) = Olig (slice start end dna) start end

--excluding end
slice :: Int -> Int -> [a] -> [a]
slice start end xs | start < 0 || end < 0 || start > end = error "incorrect coordinates"
                   | otherwise = take (end - start) (drop start xs)
