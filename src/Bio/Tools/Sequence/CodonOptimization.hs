module Bio.Tools.Sequence.CodonOptimization (
    optimizeCodonForAA
   ,optimizeCodonForDNA
   ,score
   ,CodonOptimizationConfig(..)
) where

import           Bio.Tools.Sequence.CodonOptimization.Algo  (optimizeCodonForAA,
                                                             optimizeCodonForDNA,
                                                             score)
import           Bio.Tools.Sequence.CodonOptimization.Types (CodonOptimizationConfig (..))

