module Bio.Tools.Sequence.OligoDesigner.Utils.MutationUtils (
  weightedRandom
 ,randomCodon
 ,oneMutation
 ,mutate
 ,mutateSlice
) where

import Bio.NucleicAcid.Nucleotide.Type (DNA)
import System.Random (StdGen, randomR)
import Control.Monad.State.Lazy (State)
import Bio.Tools.Sequence.OligoDesigner.Types (Codon, Weight)
import Bio.Tools.Sequence.CodonOptimization.Types (Organism)
import           Data.Map                                       as Map (lookup)
import Data.Maybe (fromMaybe)
import Bio.Protein.AminoAcid (AA(..))
import Bio.Tools.Sequence.CodonOptimization.Constants (ak2Codon, codonFrequencies, codon2ak)
import Data.List (sortOn, nub)
import Control.Monad.State (get, put, gets)
import System.Random.Shuffle (shuffle')
import Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils (slice)

oneMutation :: Organism -> Codon -> State StdGen [DNA]
oneMutation organism codon = do
    let aa = fromMaybe (error ("cannot find aa for codon " ++ show codon)) (Map.lookup codon codon2ak)
    newCodon <- randomCodon organism aa
    if newCodon == codon && aa /= TRP && aa /= MET
        then oneMutation organism codon
        else return newCodon

randomCodon :: Organism -> AA -> State StdGen Codon
randomCodon organism aa = do
    let codons = fromMaybe [] (Map.lookup aa ak2Codon)
    let codonToWeight = map (\codon -> (codon, codonWeight codon)) codons
    weightedRandom codonToWeight
  where
    codonWeight :: Codon -> Weight
    codonWeight codon = fromMaybe 0 (Map.lookup codon (codonFrequencies organism))

weightedRandom :: [(a, Weight)] -> State StdGen a
weightedRandom []    = error "cannot get random for empty array" --FIXME: normal exception instead
weightedRandom items = do
    gen <- get
    let (random, newGen) = randomR (0, sum $ map snd items) gen
    put newGen
    let value = getWeightedItem random 0 (sortOn snd items)
    return value
  where
    getWeightedItem :: Double -> Double -> [(a, Weight)] -> a
    getWeightedItem _ _ []       = error "some strange situation"
    getWeightedItem _ _ [(x, _)] = x
    getWeightedItem randomValue acc ((x, w) : xs) | acc + w >= randomValue = x
                                                  | otherwise = getWeightedItem randomValue (acc + w) xs

mutate :: Organism -> [DNA] -> (Int, Int) -> State StdGen [[DNA]]
mutate organism dna interval@(start, end) | validateInterval interval (length dna) = error ("invalid interval for mutation: " ++ show interval)
                                          | otherwise = do
    let sliceIndex = (start - 1) * 3
    let sliceEndIndex = (end - 1) * 3 + 3
    let begin = take sliceIndex dna
    let mutated = slice sliceIndex sliceEndIndex dna
    let final = drop sliceEndIndex dna
    variants <- mutateSlice organism mutated
    return $ map (\var -> begin ++ var ++ final) variants

validateInterval :: (Int, Int) -> Int -> Bool
validateInterval (start, end) len = start > end || start < 0 || end < 0 || start * 3 > len || end * 3 > len

mutateSlice :: Organism -> [DNA] -> State StdGen [[DNA]]
mutateSlice organism mutated = mutateEachCodon 0 [mutated]
  where
    mutateEachCodon :: Int -> [[DNA]] -> State StdGen [[DNA]]
    mutateEachCodon index acc
        | index * 3 >= length mutated = gets (nub . shuffle' acc (length acc))
        | otherwise = do
            let codonIndex = index * 3
            let codonEndIndex = codonIndex + 3
            let codon = take 3 (drop codonIndex mutated)
            newCodon <- oneMutation organism codon
            let variant = take codonIndex mutated ++ newCodon ++ drop codonEndIndex mutated
            mutateEachCodon (index + 1) (acc ++ [variant])