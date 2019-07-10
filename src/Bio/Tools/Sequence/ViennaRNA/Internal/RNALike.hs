module Bio.Tools.Sequence.ViennaRNA.Internal.RNALike
  ( RNALike (..)
  ) where

import           Bio.NucleicAcid.Nucleotide (DNA, RNA)
import qualified Bio.NucleicAcid.Nucleotide as N (toRNA)

-- | Class that describes objects that can be converted to RNA.
--
class RNALike a where
  toRNA :: a -> RNA

instance RNALike DNA where
  toRNA = N.toRNA

instance RNALike RNA where
  toRNA = id
