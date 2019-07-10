module Bio.Tools.Sequence.Primers.Properties
 ( calculateTargetToFold
 ) where

import           Bio.NucleicAcid.Nucleotide           (Complementary (..))
import           Bio.Tools.Sequence.Primers.Constants (annealingTemp)
import           Bio.Tools.Sequence.Primers.Types     (Primer)
import           Bio.Tools.Sequence.ViennaRNA.Cofold  (cofold)
import           Bio.Tools.Sequence.ViennaRNA.Fold    (fold)

-- | For given primer calculates relation of this primer's energy of interaction
-- with target to its energy of forming secondary structure.
--
calculateTargetToFold :: Primer -> Float
calculateTargetToFold primer = tgtEnergy / foldingEnergy
  where
    tgtEnergy     = abs $ fst $ cofold annealingTemp (reverse primer, fmap cNA primer)
    foldingEnergy = abs $ fst $ fold annealingTemp primer
