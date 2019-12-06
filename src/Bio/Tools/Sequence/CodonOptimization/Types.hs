module Bio.Tools.Sequence.CodonOptimization.Types
    ( CodonScoreConfig(..)
    , Organism(..)
    , standardTemperature
    , defaultForbiddenRegexp
    ) where

import qualified Data.ByteString.Lazy as BSL (ByteString)
import           Data.Default                (Default (..))
import           GHC.Generics                (Generic)

standardTemperature :: Double
standardTemperature = 37

data Organism = CHO | EColi | Human
    deriving (Eq, Show, Generic)

-- | all parameters for codon optimization
data CodonScoreConfig =
    CodonScoreConfig
        { organism           :: Organism
        , initLen            :: Int                -- ^ number of first ak from initial sequence, which will optimised without scoring function
        , windowLen          :: Int                -- ^ length of variation window
        , codonUsageWeight   :: Double             -- ^ Codon usage weight
        , gcWeight           :: Double             -- ^ GC-content weight
        , gcFactor           :: Double             -- ^ GC_score in the power of F_gc is used
        , gcWindow           :: Int                -- ^ length of the window for GC-score calculation (bp)
        , rnaFoldingWeight   :: Float              -- ^ Weight of the RNA folding score
        , rnaFoldingFactor   :: Float              -- ^ RNA folding score in the power of F_rnaf is used
        , rnaFoldingWindow   :: Int                -- ^ length of the window for RNA folding score calculation (bp)
        , forbiddenDNAWeight :: Double             -- ^ forbidden DNA motifs score weight
        , gcContentDesired   :: Int                -- ^ desired gc content in percents
        , forbiddenSequence  :: [BSL.ByteString]   -- ^ list of forbidden patterns
        }
    deriving (Eq, Show, Generic)

instance Default CodonScoreConfig where
  def = CodonScoreConfig CHO 3 1 1 0.5 1.4 40 0.001 2.6 100 1 43 defaultForbiddenRegexp

defaultForbiddenRegexp :: [BSL.ByteString]
defaultForbiddenRegexp =
    [ "ATTTA"
    , "ATACTCCCCC"
    , "CGATCG" -- PvuI
    , "GGGGACTTTGCACTGGAACTTACAACACCCCAGCAAGGACGCG"
    , "CCGGCGGGT"
    , "TTTATAATTTCTTCTTCCAGAA"
    , "CCGTGCTGGCGTCTG"
    , "AATAAA.{10,30}CA{30,}(TCTG|TG.CT)"
    , "(GGG|CCA)CGCCTATAAA(((C|T)(C|T)A.(T|A)(C|T)(C|T))|(TCA(G|T)T(T|C)))(A|G)G(A|T)(C|T)(G|A|C)"
    , "CAGG"
    , "(A|C)AGGT(A|G)AGT"
    , "AATAAA"
    , "GCC(A|G)CCATGG"
    ]
