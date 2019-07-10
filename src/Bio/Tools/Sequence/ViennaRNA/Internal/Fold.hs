module Bio.Tools.Sequence.ViennaRNA.Internal.Fold
    ( fold
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

foreign import ccall "vrna_fold_temperature" vrna_fold_temperature :: CString -> CString -> CDouble -> CFloat

vRnaFoldString :: Double -> String -> (Float, String)
vRnaFoldString temperature rnaSequence = unsafePerformIO $ do
    cRnaSequence <- newCString rnaSequence
    cStructure <- newCString rnaSequence -- use the same rnaSequence just to get string of the same size
    let energy = realToFrac $ vrna_fold_temperature cRnaSequence cStructure (realToFrac temperature)
    structureRes <- energy `seq` peekCString cStructure
    free cStructure
    free cRnaSequence
    return (energy, structureRes)

-- | Calculates folding energy of given 'RNALike' strand. Also returns secondary
-- structure of that strand in dot-plot form.
--
fold :: RNALike a => Double -> [a] -> (Float, String)
fold temperature = vRnaFoldString temperature . (symbol . toRNA <$>)
