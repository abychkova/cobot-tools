module Bio.Tools.Sequence.OligoDesigner.Optimizer.RNACofoldOptimizer
    ( rnaOptimize
    , maxPairMutationIndexes
    , minPairMutationIndexes
    , mutationIndexes
    ) where

import           Bio.NucleicAcid.Nucleotide.Type            (DNA)
import           Bio.Tools.Sequence.CodonOptimization       (CodonOptimizationConfig (..))
import           Bio.Tools.Sequence.CodonOptimization.Types (Organism)
import           Bio.Tools.Sequence.OligoDesigner.Scorer    (rnaScore, rnaMatrixScore)
import           Bio.Tools.Sequence.OligoDesigner.Types     (MatrixCell (..),
                                                             Olig (..),
                                                             OligSet (..),
                                                             OligsDesignerConfig (..),
                                                             OligsDesignerConfig (..), standardTemperature, OligLight(..))
import           Bio.Tools.Sequence.OligoDesigner.Utils     (assemble,
                                                             buildOligSet,
                                                             oneMutation, slice, mutateSlice, getAAIndex, mutate, compareBySecond)
import           Control.Monad.State                        (State)
import           Data.Foldable                              (minimumBy)
import           Data.List                                  (intersect,
                                                             maximumBy, nub, findIndex)
import           Data.Matrix                                (Matrix, ncols,
                                                             nrows, (!), matrix)
import           System.Random                              (StdGen)
import Bio.Tools.Sequence.ViennaRNA.Internal.Cofold (cofold)
import Text.Regex.TDFA (Regex)
import Debug.Trace
import Bio.Tools.Sequence.OligoDesigner.RNAMatrixBuilder (rnaMatrix, rebuildMatrix)
import Bio.Tools.Sequence.OligoDesigner.Prettifier (prettyDNA, prettyMatrixCell, prettyOligSet)
import Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer (filterForbidden)

rnaOptimize :: OligsDesignerConfig -> [Regex] -> OligSet -> State StdGen OligSet
rnaOptimize conf regexes oligs@(OligSet _ _ splitting) = do
    let organismType = organism $ codonOptimizationConfig conf
    let mtx = rnaMatrix oligs
    let indexesToMutate = mutationIndexes mtx
    let dna = assemble oligs
    sequenceVariants <- concat <$> mapM (mutate organismType dna) indexesToMutate
    let filtered =  filterForbidden regexes sequenceVariants
    let oligsVariants = map (buildOligSet splitting) filtered
    let oligs2score = map (scoreOligs mtx) oligsVariants
    return $ fst $ maximumBy compareBySecond oligs2score
  where
    scoreOligs :: Matrix MatrixCell -> OligSet -> (OligSet, Float)
    scoreOligs mtx oligs = (oligs, score)
      where
        score = rnaMatrixScore $ rebuildMatrix mtx oligs

mutationIndexes :: Matrix MatrixCell -> [(Int, Int)]
mutationIndexes oligsMatrix = nub (minPairMutationIndexes minPair ++ maxPairMutationIndexes maxPair)
  where
    rowsCnt = nrows oligsMatrix
    colsCnt = ncols oligsMatrix

    minPair = minimumBy compareByRna [oligsMatrix ! (x , x + 1) | x <- [1 .. rowsCnt - 1]]
    maxPair = maximumBy compareByRna [oligsMatrix ! (x, y) | x <- [1 .. rowsCnt], y <- [1 .. colsCnt], x <= y && abs (x - y) /= 1]

    compareByRna :: MatrixCell -> MatrixCell -> Ordering
    compareByRna (MatrixCell _ _ rna1) (MatrixCell _ _ rna2) = compare (abs rna1) (abs rna2)

minPairMutationIndexes :: MatrixCell -> [(Int, Int)]
minPairMutationIndexes (MatrixCell (OligLight _ (Olig _ start1 end1)) (OligLight _ (Olig _ start2 end2)) _) = indexes
  where
    intersection = [start1 .. end1 - 1] `intersect` [start2 .. end2 - 1]
    indexes =
        if null intersection
            then []
            else [(getAAIndex (head intersection), getAAIndex (last intersection))]

maxPairMutationIndexes :: MatrixCell -> [(Int, Int)]
maxPairMutationIndexes (MatrixCell (OligLight _ (Olig _ start1 end1)) (OligLight _ (Olig _ start2 end2)) _) =
    [(getAAIndex start1, getAAIndex $ end1 - 1), (getAAIndex start2, getAAIndex $ end2 - 1)]