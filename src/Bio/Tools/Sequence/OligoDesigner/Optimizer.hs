module Bio.Tools.Sequence.OligoDesigner.Optimizer
    ( minMaxOptimize
    , maxPairMutationIndexes
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
                                                             oneMutation, slice)
import           Control.Monad.State                        (State)
import           Data.Foldable                              (minimumBy)
import           Data.List                                  (intersect,
                                                             maximumBy)
import           Data.Matrix                                (Matrix, ncols,
                                                             nrows, (!))
import           System.Random                              (StdGen)

--TODO: test me (min-max optimizer)
minMaxOptimize :: OligoDesignerConfig -> OligSet -> State StdGen OligSet
minMaxOptimize conf@(OligoDesignerConfig codonConf _ _) oligs@(OligSet _ _ splitting) = do
    let dna = assemble oligs
    let indexesToMutate = mutationIndexes $ rnaMatrix oligs
    let mutateFunction = mutate codonConf dna
    varSequences <- concat <$> mapM mutateFunction indexesToMutate
    let variants = map (buildOligSet splitting) varSequences
    return $ maximumBy scoreCmp variants
  where
    scoreCmp :: OligSet -> OligSet -> Ordering
    scoreCmp oligs1 oligs2 = compare score1 score2
      where
        score1 = commonScore conf oligs1
        score2 = commonScore conf oligs2

--TODO: test me
mutationIndexes :: Matrix MatrixCell -> [(Int, Int)]
mutationIndexes oligsMatrix = minPairMutationIndexes minPair ++ maxPairMutationIndexes maxPair
  where
    rowsCnt = nrows oligsMatrix
    colsCnt = ncols oligsMatrix
    minPair = minimumBy orderByRna [oligsMatrix ! (x, y) | x <- [1 .. rowsCnt], y <- [1 .. colsCnt], abs (x - y) == 1]
    maxPair = maximumBy orderByRna [oligsMatrix ! (x, y) | x <- [1 .. rowsCnt], y <- [1 .. colsCnt], abs (x - y) /= 1]
    
    orderByRna :: MatrixCell -> MatrixCell -> Ordering
    orderByRna (MatrixCell _ _ rna1) (MatrixCell _ _ rna2) = compare rna1 rna2

--TODO: test me
minPairMutationIndexes :: MatrixCell -> [(Int, Int)]
minPairMutationIndexes (MatrixCell (Olig _ start1 end1) (Olig _ start2 end2) _) = indexes
  where
    intersection = [start1 .. end1] `intersect` [start2 .. end2]
    indexes =
        if null intersection
            then []
            else [(div (head intersection) 3, div (last intersection) 3)]

--TODO: test me
maxPairMutationIndexes :: MatrixCell -> [(Int, Int)]
maxPairMutationIndexes (MatrixCell (Olig _ start1 end1) (Olig _ start2 end2) _) =
    [(div start1 3, div end1 3), (div start2 3, div end2 3)]

--TODO: test me
mutate :: CodonOptimizationConfig -> [DNA] -> (Int, Int) -> State StdGen [[DNA]]
mutate conf dna (start, end) = do
    let sliceIndex = start * 3
    let sliceEndIndex = end * 3 + 3
    let begin = take sliceIndex dna
    let mutual = slice sliceIndex sliceEndIndex dna
    let final = drop sliceEndIndex dna
    variants <- mutateSlice (organism conf) mutual
    return $ map (\var -> begin ++ var ++ final) variants

--TODO: test me
mutateSlice :: Organism -> [DNA] -> State StdGen [[DNA]]
mutateSlice organism dna = mutateEachCodon dna 0 [dna]
  where
    mutateEachCodon :: [DNA] -> Int -> [[DNA]] -> State StdGen [[DNA]]
    mutateEachCodon mutual index acc
        | index * 3 >= length mutual = return acc
        | otherwise = do
            let codonIndex = index * 3
            let codonEndIndex = codonIndex + 3
            let codon = take 3 (drop codonIndex mutual)
            newCodon <- oneMutation organism codon
            let variant = take codonIndex mutual ++ newCodon ++ drop codonEndIndex mutual
            mutateEachCodon mutual (index + 1) (acc ++ [variant])
