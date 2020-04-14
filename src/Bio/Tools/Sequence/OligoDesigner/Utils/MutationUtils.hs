module Bio.Tools.Sequence.OligoDesigner.Utils.MutationUtils
    ( weightedRandom
    , randomCodon
    , oneMutation
    , mutate
    , mutateSlice
    ) where

import Control.Monad.Except           (Except, throwError)
import Control.Monad.Trans            (lift)
import Control.Monad.Trans.State.Lazy (StateT, get, gets, put)
import Data.List                      (nub, sortOn)
import Data.Map                       as Map (lookup)
import Data.Maybe                     (fromMaybe)
import System.Random                  (StdGen, randomR)
import System.Random.Shuffle          (shuffle')

import Bio.NucleicAcid.Nucleotide.Type (DNA)
import Bio.Protein.AminoAcid           (AA (..))

import Bio.Tools.Sequence.CodonOptimization.Constants     (ak2Codon, codon2ak, codonFrequencies)
import Bio.Tools.Sequence.CodonOptimization.Types         (Organism)
import Bio.Tools.Sequence.OligoDesigner.Types             (Codon, Weight)
import Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils (slice)


mutate :: Organism -> [DNA] -> (Int, Int) -> StateT StdGen (Except String) [[DNA]]
mutate organism dna interval@(start, end)
  | validateInterval interval (length dna) = throwError ("invalid interval for mutation: " ++ show interval)
  | otherwise = do
    let sliceIndex = (start - 1) * 3
    let sliceEndIndex = (end - 1) * 3 + 3
    let begin = take sliceIndex dna
    mutated <- lift $ slice sliceIndex sliceEndIndex dna
    let final = drop sliceEndIndex dna
    variants <- mutateSlice organism mutated
    return $ map (\var -> begin ++ var ++ final) variants

validateInterval :: (Int, Int) -> Int -> Bool
validateInterval (start, end) len = start > end || start < 0 || end < 0 || start * 3 > len || end * 3 > len

mutateSlice :: Organism -> [DNA] -> StateT StdGen (Except String) [[DNA]]
mutateSlice organism mutated = mutateEachCodon 0 [mutated]
  where
    mutateEachCodon :: Int -> [[DNA]] -> StateT StdGen (Except String) [[DNA]]
    mutateEachCodon index acc
        | index * 3 >= length mutated = gets (nub . shuffle' acc (length acc))
        | otherwise = do
            let codonIndex = index * 3
            let codonEndIndex = codonIndex + 3
            let codon = take 3 (drop codonIndex mutated)
            newCodon <- oneMutation organism codon
            let variant = take codonIndex mutated ++ newCodon ++ drop codonEndIndex mutated
            mutateEachCodon (index + 1) (acc ++ [variant])

oneMutation :: Organism -> Codon -> StateT StdGen (Except String) Codon
oneMutation organism codon =
    case Map.lookup codon codon2ak of
        Nothing -> throwError ("cannot find aa for codon " ++ show codon)
        Just aa -> do
            newCodon <- randomCodon organism aa
            if newCodon == codon && aa /= TRP && aa /= MET
                then oneMutation organism codon
                else return newCodon

randomCodon :: Organism -> AA -> StateT StdGen (Except String) Codon
randomCodon organism aa = do
    let codons = fromMaybe [] (Map.lookup aa ak2Codon)
    let codonToWeight = map (\codon -> (codon, codonWeight codon)) codons
    weightedRandom codonToWeight
  where
    codonWeight :: Codon -> Weight
    codonWeight codon = fromMaybe 0 (Map.lookup codon (codonFrequencies organism))

weightedRandom :: [(a, Weight)] -> StateT StdGen (Except String) a
weightedRandom []    = throwError "cannot get random for empty array"
weightedRandom items = do
    gen <- get
    let (random, newGen) = randomR (0, sum $ map snd items) gen
    put newGen
    lift $ getWeightedItem random 0 (sortOn snd items)
  where
    getWeightedItem :: Double -> Double -> [(a, Weight)] -> Except String a
    getWeightedItem _ _ []       = throwError "some strange situation"
    getWeightedItem _ _ [(x, _)] = return x
    getWeightedItem randomValue acc ((x, w) : xs) 
         | acc + w >= randomValue = return x
         | otherwise = getWeightedItem randomValue (acc + w) xs
