module Bio.Tools.Sequence.ViennaRNA.Internal.Cofold
    ( cofold
    ) where

import           Bio.NucleicAcid.Nucleotide                    (symbol)
import           Bio.Tools.Sequence.ViennaRNA.Internal.RNALike (RNALike (..))
import           Foreign.C.String                              (CString,
                                                                newCString,
                                                                peekCString)
import           Foreign.C.Types                               (CDouble (..),
                                                                CFloat (..))
import           Foreign.Marshal.Alloc                         (free)
import           System.IO.Unsafe                              (unsafePerformIO)

foreign import ccall "vrna_cofold_temperature" vrna_cofold_temperature :: CString -> CString -> CDouble -> CFloat

-- TODO make it for list of pairs to allocate CString of only one size
-- TODO (need to handle different sizes of oligs, try to take just max size (because CString ends with \0)
vRnaCofoldString :: Double -> (String, String) -> (Float, String)
vRnaCofoldString temperature (rnaSequence1, rnaSequence2) = unsafePerformIO $ do
    let inputSeq = rnaSequence1 ++ "&" ++ rnaSequence2
    cRnaSequence <- newCString' inputSeq
    let energy = realToFrac $ vrna_cofold_temperature cRnaSequence cRnaSequence (realToFrac temperature)
    structureRes <- evalResult energy cRnaSequence
    free' cRnaSequence
    return (energy, structureRes)
    where
      newCString' = newCString
      free' = free
      evalResult energy cRnaSequence = energy `seq` peekCString cRnaSequence

-- | Calculates cofolding energy and interaction dot-plot between two 'RNALike' strands.
--
cofold :: RNALike a => Double -> ([a], [a]) -> (Float, String)
cofold temperature = vRnaCofoldString temperature . toStringsPair
  where
    toStringsPair (nucs1, nucs2) = (symbol . toRNA <$> nucs1, symbol . toRNA <$> nucs2)
