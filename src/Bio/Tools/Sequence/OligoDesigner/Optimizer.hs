module Bio.Tools.Sequence.OligoDesigner.Optimizer
    ( minMaxOptimize
    , maxPairMutationIndexes
    , minPairMutationIndexes
    , mutationIndexes
    , mutateSlice
    , mutate
    ) where

import           Bio.NucleicAcid.Nucleotide.Type            (DNA)
import           Bio.Tools.Sequence.CodonOptimization       (CodonOptimizationConfig (..))
import           Bio.Tools.Sequence.CodonOptimization.Types (Organism)
import           Bio.Tools.Sequence.OligoDesigner.Scorer    (commonScore,
                                                             rnaMatrix)
import           Bio.Tools.Sequence.OligoDesigner.Types     (MatrixCell (..),
                                                             Olig (..),
                                                             OligSet (..),
                                                             OligoDesignerConfig (..),
                                                             OligoDesignerConfig (..))
import           Bio.Tools.Sequence.OligoDesigner.Utils     (assemble,
                                                             buildOligSet,
                                                             oneMutation, slice, prettyDNA)
import           Control.Monad.State                        (State)
import           Data.Foldable                              (minimumBy)
import           Data.List                                  (intersect,
                                                             maximumBy, nub)
import           Data.Matrix                                (Matrix, ncols,
                                                             nrows, (!), prettyMatrix)
import           System.Random                              (StdGen)
import Debug.Trace (trace)

minMaxOptimize :: OligoDesignerConfig -> OligSet -> State StdGen OligSet
minMaxOptimize conf@(OligoDesignerConfig codonConf _ _) oligs@(OligSet _ _ splitting) = do
    let dna = assemble oligs
    let indexesToMutate = mutationIndexes $ rnaMatrix oligs
    let mutateFunction = mutate codonConf dna
    varSequences <- concat <$> mapM mutateFunction indexesToMutate
    let variants = map (buildOligSet splitting) varSequences
    let max = maximumBy scoreCmp variants
    return max
  where
    scoreCmp :: OligSet -> OligSet -> Ordering
    scoreCmp oligs1 oligs2 = compare score1 score2
      where
        score1 = commonScore conf oligs1
        score2 = commonScore conf oligs2

mutationIndexes :: Matrix MatrixCell -> [(Int, Int)]
mutationIndexes oligsMatrix = nub (minPairMutationIndexes minPair ++ maxPairMutationIndexes maxPair)
  where
    rowsCnt = nrows oligsMatrix
    colsCnt = ncols oligsMatrix
    minPair = minimumBy compareByRna [oligsMatrix ! (x, y) | x <- [1 .. rowsCnt], y <- [1 .. colsCnt], abs (x - y) == 1]
    maxPair = maximumBy compareByRna [oligsMatrix ! (x, y) | x <- [1 .. rowsCnt], y <- [1 .. colsCnt], abs (x - y) /= 1]

    compareByRna :: MatrixCell -> MatrixCell -> Ordering
    compareByRna (MatrixCell _ _ rna1) (MatrixCell _ _ rna2) = compare rna1 rna2

minPairMutationIndexes :: MatrixCell -> [(Int, Int)]
minPairMutationIndexes (MatrixCell (Olig _ start1 end1) (Olig _ start2 end2) _) = indexes
  where
    intersection = [start1 .. end1 - 1] `intersect` [start2 .. end2 - 1]
    indexes =
        if null intersection
            then []
            else [(getAANumber (head intersection), getAANumber (last intersection))]

maxPairMutationIndexes :: MatrixCell -> [(Int, Int)]
maxPairMutationIndexes (MatrixCell (Olig _ start1 end1) (Olig _ start2 end2) _) =
    [(getAANumber start1, getAANumber $ end1 - 1), (getAANumber start2, getAANumber $ end2 - 1)]

getAANumber :: Int -> Int
getAANumber coordinate = ceiling (realToFrac (coordinate + 1) / 3)

mutate :: CodonOptimizationConfig -> [DNA] -> (Int, Int) -> State StdGen [[DNA]]
mutate conf dna interval@(start, end) | validateInterval interval (length dna) = error ("invalid interval for mutation: " ++ show interval)
                                      | otherwise = do
    let sliceIndex = (start - 1) * 3
    let sliceEndIndex = (end - 1) * 3 + 3
    let begin = take sliceIndex dna
    let mutated = slice sliceIndex sliceEndIndex dna
    let final = drop sliceEndIndex dna
    variants <- mutateSlice (organism conf) mutated
    return $ map (\var -> begin ++ var ++ final) variants
    
mutateSlice :: Organism -> [DNA] -> State StdGen [[DNA]]
mutateSlice organism dna = mutateEachCodon dna 0 [dna]
  where
    mutateEachCodon :: [DNA] -> Int -> [[DNA]] -> State StdGen [[DNA]]
    mutateEachCodon mutated index acc
        | index * 3 >= length mutated = return $ nub acc
        | otherwise = do
            let codonIndex = index * 3
            let codonEndIndex = codonIndex + 3
            let codon = take 3 (drop codonIndex mutated)
            newCodon <- oneMutation organism codon
            let variant = take codonIndex mutated ++ newCodon ++ drop codonEndIndex mutated
            mutateEachCodon mutated (index + 1) (acc ++ [variant])

validateInterval :: (Int, Int) -> Int -> Bool
validateInterval (start, end) len = start > end || start < 0 || end < 0 || start * 3 > len || end * 3 > len