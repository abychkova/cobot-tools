module Bio.Tools.Sequence.OligoDesigner.Utils
 (assemble
 ,weightedRandom
 ,randomCodon
 ,buildOligSet
 ,slice
 ,translate
 ,oneMutation
 ,mutate
 ,mutateSlice
 ,getAAIndex
 ,mixOligs
 ,compareBySecond
 ) where

import           Bio.NucleicAcid.Nucleotide.Type                (DNA (..), cNA)
import           Bio.Protein.AminoAcid.Type                     (AA (..))
import           Bio.Tools.Sequence.CodonOptimization.Constants (ak2Codon,
                                                                 codon2ak,
                                                                 codonFrequencies)
import           Bio.Tools.Sequence.CodonOptimization.Types     (Organism)
import           Bio.Tools.Sequence.OligoDesigner.Types         (Codon,
                                                                 Olig (..),
                                                                 OligBounds,
                                                                 OligSet (..),
                                                                 OligSplitting (..), OligsDesignerConfig)
import           Control.Monad.State                            (State, get,
                                                                 put)
import           Data.List                                      (sortOn, nub)
import           Data.Map                                       as Map (lookup)
import           Data.Maybe                                     (fromMaybe)
import           System.Random                                  (StdGen,
                                                                 randomR)
import Debug.Trace 
import Bio.Tools.Sequence.CodonOptimization (CodonOptimizationConfig(..))
import Text.Regex.TDFA (Regex, makeRegex, match)
import Control.Monad.State.Lazy (gets)
import System.Random.Shuffle (shuffle')
        
mixOligs :: OligSet -> [Olig]
mixOligs (OligSet forward reversed _) = mix forward reversed
  where
    mix :: [a] -> [a] -> [a]
    mix (x:xs) (y:ys) = x : y : mix xs ys
    mix x []          = x
    mix [] y          = y
    
assemble :: OligSet -> [DNA]
assemble (OligSet fwd rvd _) = construct fwd rvd 0 [] where

    construct :: [Olig] -> [Olig] -> Int -> [DNA] -> [DNA]
    construct _ [] _ acc = acc
    construct [] _ _ acc = acc
    construct (Olig seq1 startLeft endLeft : xs) (Olig seq2 startRight endRight : ys) prevEnd acc = construct xs ys endRight res
      where
        partFormOlig1 = drop (prevEnd - startLeft) seq1
        partFormOlig2 = map cNA (drop (endLeft - startRight) (reverse seq2))
        res = acc ++ partFormOlig1 ++ partFormOlig2

--TODO: correct exception instead
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
    strand5' = map (buildOlig id sequ) (strand5 splitting)
    strand3' = map (buildOlig reverse (map cNA sequ)) (strand3 splitting)

    buildOlig :: ([DNA] -> [DNA]) -> [DNA] -> OligBounds -> Olig
    buildOlig fun dna (start, end) = Olig sliceDNA start end
      where
        sliceDNA = fun $ slice start end dna

--excluding end
slice :: Int -> Int -> [a] -> [a]
slice start end xs | start < 0 || end < 0 || start > end = error "incorrect coordinates"
                   | otherwise = take (end - start) (drop start xs)

translate :: [DNA] -> [DNA]
translate = map cNA

--FIXME: что если вернется пустой список вариантов, потому что во всех есть запрещенки? что сказать пользователю?
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

--TODO: test me  
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
            
--TODO: test me      
getAAIndex :: Int -> Int
getAAIndex coordinate = ceiling (realToFrac (coordinate + 1) / 3)

compareBySecond :: Ord b => (a, b) -> (a, b) -> Ordering
compareBySecond p1 p2 = compare (snd p1) (snd p2)