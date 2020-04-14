module Bio.Tools.Sequence.OligoDesigner.Optimizer.RNACofoldOptimizer
    ( rnaOptimize
    , maxPairMutationIndexes
    , minPairMutationIndexes
    , mutationIndexes
    ) where

import           Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer         (filterForbidden)
import           Bio.Tools.Sequence.OligoDesigner.Scorer                 (rnaMatrixScore)
import           Bio.Tools.Sequence.OligoDesigner.Types                  (MatrixCell (..),
                                                                          Olig (..),
                                                                          OligLight (..),
                                                                          OligSet (..),
                                                                          OligsDesignerInnerConfig (..))
import           Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils      (assemble,
                                                                          buildOligSet,
                                                                          getAAIndex,
                                                                          orderByScore)
import           Bio.Tools.Sequence.OligoDesigner.Utils.MutationUtils    (mutate)
import           Bio.Tools.Sequence.OligoDesigner.Utils.RNAMatrixBuilder (rebuildMatrix,
                                                                          rnaMatrix)
import           Control.Monad.Except                                    (Except)
import           Control.Monad.Trans                                     (lift)
import           Control.Monad.Trans.State.Lazy                          (StateT)
import           Data.Foldable                                           (minimumBy)
import           Data.List                                               (intersect,
                                                                          maximumBy,
                                                                          nub)
import           Data.Matrix                                             (Matrix,
                                                                          ncols,
                                                                          nrows,
                                                                          (!))
import           System.Random                                           (StdGen)

-- | 'rnaOptimize' function does optimization for oligs rna matrix.
-- The main idea is to reduce the difference between neighboring oligs with minimum rna energy
-- and not neighboring oligs with maximum rna energy
rnaOptimize :: OligsDesignerInnerConfig              -- ^ configuration data (used 'organism' and 'regexes')
            -> OligSet                               -- ^ oligs set
            -> StateT StdGen (Except String) OligSet -- ^ result of optimization is optimized oligs set with the best score or error string
rnaOptimize (OligsDesignerInnerConfig organism _ regexes _ _) oligs@(OligSet _ _ splitting) = do
    let mtx = rnaMatrix oligs
    let dna = assemble oligs
    indexesToMutate <- lift $ mutationIndexes mtx
    sequenceVariants <- concat <$> mapM (mutate organism dna) indexesToMutate
    let filtered = filterForbidden regexes sequenceVariants
    oligsVariants <- lift $ mapM (buildOligSet splitting) filtered
    let (_, maxOligs) = orderByScore oligsVariants (rnaMatrixScore . rebuildMatrix mtx)
    return maxOligs

mutationIndexes :: Matrix MatrixCell -> Except String [(Int, Int)]
mutationIndexes oligsMatrix = do
    let rowsCnt = nrows oligsMatrix
    let colsCnt = ncols oligsMatrix
    let minPair = minimumBy compareByRna [oligsMatrix ! (x, x + 1) | x <- [1 .. rowsCnt - 1]]
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
    endAA <- getAAIndex end
    return (startAA, endAA)