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
                                                             OligsDesignerInnerConfig (..), standardTemperature, OligLight(..))
import           Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils     (assemble,
                                                             buildOligSet,
                                                             slice, getAAIndex, compareBySecond, orderByScore)
import Bio.Tools.Sequence.OligoDesigner.Utils.MutationUtils (oneMutation, mutateSlice, mutate)
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
import Bio.Tools.Sequence.OligoDesigner.Utils.RNAMatrixBuilder (rnaMatrix, rebuildMatrix)
import Bio.Tools.Sequence.OligoDesigner.Utils.Prettifier (prettyDNA, prettyMatrixCell, prettyOligSet)
import Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer (filterForbidden)
import Control.Monad.Trans.State.Lazy (StateT)
import Control.Monad.Except (Except)
import Control.Monad.Trans (lift)

rnaOptimize :: OligsDesignerInnerConfig -> OligSet -> StateT StdGen (Except String) OligSet
rnaOptimize (OligsDesignerInnerConfig organism _ regexes _ _) oligs@(OligSet _ _ splitting) = do
    let mtx = rnaMatrix oligs
    indexesToMutate <- lift $ mutationIndexes mtx
    let dna = assemble oligs
    sequenceVariants <- concat <$> mapM (mutate organism dna) indexesToMutate
    let filtered = filterForbidden regexes sequenceVariants
    oligsVariants <- lift $ mapM (buildOligSet splitting) filtered
    let (_, max) = orderByScore oligsVariants (rnaMatrixScore . rebuildMatrix mtx)
    return max

mutationIndexes :: Matrix MatrixCell -> Except String [(Int, Int)]
mutationIndexes oligsMatrix = do
    let rowsCnt = nrows oligsMatrix
    let colsCnt = ncols oligsMatrix

    let minPair = minimumBy compareByRna [oligsMatrix ! (x , x + 1) | x <- [1 .. rowsCnt - 1]]
    let maxPair = maximumBy compareByRna [oligsMatrix ! (x, y) | x <- [1 .. rowsCnt], y <- [1 .. colsCnt], x <= y && abs (x - y) /= 1]

    minPairIndexes <- minPairMutationIndexes minPair
    maxPairIndexes <- maxPairMutationIndexes maxPair

    return $ nub (minPairIndexes ++ maxPairIndexes)
  where
    compareByRna :: MatrixCell -> MatrixCell -> Ordering
    compareByRna (MatrixCell _ _ rna1) (MatrixCell _ _ rna2) = compare (abs rna1) (abs rna2)

minPairMutationIndexes :: MatrixCell -> Except String [(Int, Int)]
minPairMutationIndexes (MatrixCell (OligLight _ (Olig _ start1 end1)) (OligLight _ (Olig _ start2 end2)) _) = do
    let intersection = [start1 .. end1 - 1] `intersect` [start2 .. end2 - 1]
    if null intersection
        then return []
        else (: []) <$> getAAIndexes (head intersection, last intersection)

maxPairMutationIndexes :: MatrixCell -> Except String [(Int, Int)]
maxPairMutationIndexes (MatrixCell (OligLight _ (Olig _ start1 end1)) (OligLight _ (Olig _ start2 end2)) _) = do
    indexesAA1 <- getAAIndexes (start1, end1 - 1)
    indexesAA2 <- getAAIndexes (start2, end2 - 1)
    return [indexesAA1, indexesAA2]

getAAIndexes :: (Int, Int) -> Except String (Int, Int)
getAAIndexes (start, end) = do
    startAA <- getAAIndex start
    endAA   <- getAAIndex $ end
    return (startAA, endAA)