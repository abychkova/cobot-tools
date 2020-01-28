module Bio.Tools.Sequence.OligoDesigner.RnaCofoldOptimizer
    ( minMaxOptimize
    , maxPairMutationIndexes
    , minPairMutationIndexes
    , mutationIndexes
    ) where

import           Bio.NucleicAcid.Nucleotide.Type            (DNA)
import           Bio.Tools.Sequence.CodonOptimization       (CodonOptimizationConfig (..))
import           Bio.Tools.Sequence.CodonOptimization.Types (Organism)
import           Bio.Tools.Sequence.OligoDesigner.Scorer    (rnaMatrix, rnaScore)
import           Bio.Tools.Sequence.OligoDesigner.Types     (MatrixCell (..),
                                                             Olig (..),
                                                             OligSet (..),
                                                             OligoDesignerConfig (..),
                                                             OligoDesignerConfig (..))
import           Bio.Tools.Sequence.OligoDesigner.Utils     (assemble,
                                                             buildOligSet,
                                                             oneMutation, slice, prettyDNA, mutateSlice, getAANumber, mutate)
import           Control.Monad.State                        (State)
import           Data.Foldable                              (minimumBy)
import           Data.List                                  (intersect,
                                                             maximumBy, nub, findIndex)
import           Data.Matrix                                (Matrix, ncols,
                                                             nrows, (!), prettyMatrix)
import           System.Random                              (StdGen)
import Debug.Trace (trace)

minMaxOptimize :: OligoDesignerConfig -> OligSet -> State StdGen OligSet
minMaxOptimize (OligoDesignerConfig (CodonOptimizationConfig organism _ _ _ _ _ _ _ _ _ _ _ _) _ _) oligs@(OligSet _ _ splitting) = do
    let dna = assemble oligs
    let indexesToMutate = mutationIndexes $ rnaMatrix oligs
    varSequences <- concat <$> mapM (mutate organism dna) indexesToMutate
    let variants = map (buildOligSet splitting) varSequences
    return $ maximumBy scoreCmp variants
  where
    scoreCmp :: OligSet -> OligSet -> Ordering
    scoreCmp oligs1 oligs2 = compare score1 score2
      where
        score1 = rnaScore oligs1
        score2 = rnaScore oligs2

mutationIndexes :: Matrix MatrixCell -> [(Int, Int)]
mutationIndexes oligsMatrix = nub (minPairMutationIndexes minPair ++ maxPairMutationIndexes maxPair)
  where
    rowsCnt = nrows oligsMatrix
    colsCnt = ncols oligsMatrix
    minPair = minimumBy compareByRna [oligsMatrix ! (x, y) | x <- [1 .. rowsCnt], y <- [1 .. colsCnt], abs (x - y) == 1]
    maxPair = maximumBy compareByRna [oligsMatrix ! (x, y) | x <- [1 .. rowsCnt], y <- [1 .. colsCnt], x <= y && abs (x - y) /= 1]

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