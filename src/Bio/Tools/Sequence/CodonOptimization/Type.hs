module Bio.Tools.Sequence.CodonOptimization.Type
    ( CodonConfig(..)
    , CodonScoreConfig(..)
    , ak2Codon
    , codon2ak
    , forbiddenRegexp
    , codonFrequencies
    , ak2MaxFrequCodon
    , standardTemperature
    ) where

import           Bio.NucleicAcid.Nucleotide.Type (DNA (..))
import           Bio.Protein.AminoAcid.Type      (AA (..))
import qualified Data.ByteString.Lazy            as BSL (ByteString)
import           Data.Default                    (Default (..), def)
import           Data.Map                        as Map (Map, fromList)

standardTemperature :: Double
standardTemperature = 37

-- | all parameters for codon optimization
data CodonConfig =
     CodonConfig
        { initLen   :: Int              -- ^ number of first ak from initial sequence, which will optimised without scoring function
        , windowLen :: Int              -- ^ length of variation window
        , scoreConf :: CodonScoreConfig -- ^ parameters for scoring function
        } deriving (Show)

-- | all parameters for scoring function of codon optimization
data CodonScoreConfig =
    CodonScoreConfig
        { codonUsageWeight   :: Double -- ^ Codon usage weight
        , gcWeight           :: Double -- ^ GC-content weight
        , gcFactor           :: Double -- ^ GC_score in the power of F_gc is used
        , gcWindow           :: Int    -- ^ length of the window for GC-score calculation (bp)
        , rnaFoldingWeight   :: Float  -- ^ Weight of the RNA folding score
        , rnaFoldingFactor   :: Float  -- ^ RNA folding score in the power of F_rnaf is used
        , rnaFoldingWindow   :: Int    -- ^ length of the window for RNA folding score calculation (bp)
        , forbiddenDNAWeight :: Double -- ^ forbidden DNA motifs score weight
        , gcContentDesired   :: Int    -- ^ desired gc content in percents
        } deriving (Show)

instance Default CodonScoreConfig where
  def = CodonScoreConfig 1 0.5 1.4 40 0.001 2.6 100 1 43

instance Default CodonConfig where
  def = CodonConfig 3 1 def

forbiddenRegexp :: [BSL.ByteString]
forbiddenRegexp =
    [ "ATTTA"
    , "ATACTCCCCC"
    , "CGATCG"
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

ak2Codon :: Map AA [[DNA]]
ak2Codon =
    fromList
        [ (PHE, [[DT, DT, DT], [DT, DT, DC]])
        , (TYR, [[DT, DA, DT], [DT, DA, DC]])
        , (CYS, [[DT, DG, DT], [DT, DG, DC]])
        , (TRP, [[DT, DG, DG]])
        , (LEU, [[DT, DT, DA], [DT, DT, DG], [DC, DT, DT], [DC, DT, DC], [DC, DT, DA], [DC, DT, DG]])
        , (PRO, [[DC, DC, DT], [DC, DC, DC], [DC, DC, DA], [DC, DC, DG]])
        , (HIS, [[DC, DA, DT], [DC, DA, DC]])
        , (GLN, [[DC, DA, DA], [DC, DA, DG]])
        , (ILE, [[DA, DT, DT], [DA, DT, DC], [DA, DT, DA]])
        , (MET, [[DA, DT, DG]])
        , (THR, [[DA, DC, DT], [DA, DC, DC], [DA, DC, DA], [DA, DC, DG]])
        , (ASN, [[DA, DA, DT], [DA, DA, DC]])
        , (LYS, [[DA, DA, DA], [DA, DA, DG]])
        , (SER, [[DT, DC, DT], [DT, DC, DC], [DT, DC, DA], [DT, DC, DG], [DA, DG, DT], [DA, DG, DC]])
        , (ARG, [[DC, DG, DT], [DC, DG, DC], [DC, DG, DA], [DC, DG, DG], [DA, DG, DA], [DA, DG, DG]])
        , (VAL, [[DG, DT, DT], [DG, DT, DC], [DG, DT, DA], [DG, DT, DG]])
        , (ALA, [[DG, DC, DT], [DG, DC, DC], [DG, DC, DA], [DG, DC, DG]])
        , (ASP, [[DG, DA, DT], [DG, DA, DC]])
        , (GLU, [[DG, DA, DA], [DG, DA, DG]])
        , (GLY, [[DG, DG, DT], [DG, DG, DC], [DG, DG, DA], [DG, DG, DG]])
        ]


ak2MaxFrequCodon :: Map AA ([DNA], Double)
ak2MaxFrequCodon =
    fromList
        [ (ALA, ([DG, DC, DC], 0.4))
        , (CYS, ([DT, DG, DC], 0.54))
        , (ASP, ([DG, DA, DC], 0.54))
        , (GLU, ([DG, DA, DG], 0.58))
        , (PHE, ([DT, DT, DC], 0.54))
        , (GLY, ([DG, DG, DC], 0.34))
        , (HIS, ([DC, DA, DC], 0.58))
        , (ILE, ([DA, DT, DC], 0.47))
        , (LYS, ([DA, DA, DG], 0.57))
        , (LEU, ([DC, DT, DG], 0.4))
        , (MET, ([DA, DT, DG], 1.0))
        , (ASN, ([DA, DA, DC], 0.53))
        , (PRO, ([DC, DC, DC], 0.32))
        , (GLN, ([DC, DA, DG], 0.73))
        , (ARG, ([DA, DG, DG], 0.21))
        , (SER, ([DA, DG, DC], 0.24))
        , (THR, ([DA, DC, DC], 0.36))
        , (VAL, ([DG, DT, DG], 0.46))
        , (TRP, ([DT, DG, DG], 1.0))
        , (TYR, ([DT, DA, DC], 0.56))
        ]

codon2ak :: Map [DNA] AA
codon2ak =
    fromList
        [ ([DT, DT, DT], PHE)
        , ([DT, DT, DC], PHE)
        , ([DT, DT, DA], LEU)
        , ([DT, DT, DG], LEU)
        , ([DT, DC, DT], SER)
        , ([DT, DC, DC], SER)
        , ([DT, DC, DA], SER)
        , ([DT, DC, DG], SER)
        , ([DT, DA, DT], TYR)
        , ([DT, DA, DC], TYR)
        , ([DT, DG, DT], CYS)
        , ([DT, DG, DC], CYS)
        , ([DT, DG, DG], TRP)
        , ([DC, DT, DT], LEU)
        , ([DC, DT, DC], LEU)
        , ([DC, DT, DA], LEU)
        , ([DC, DT, DG], LEU)
        , ([DC, DC, DT], PRO)
        , ([DC, DC, DC], PRO)
        , ([DC, DC, DA], PRO)
        , ([DC, DC, DG], PRO)
        , ([DC, DA, DT], HIS)
        , ([DC, DA, DC], HIS)
        , ([DC, DA, DA], GLN)
        , ([DC, DA, DG], GLN)
        , ([DC, DG, DT], ARG)
        , ([DC, DG, DC], ARG)
        , ([DC, DG, DA], ARG)
        , ([DC, DG, DG], ARG)
        , ([DA, DT, DT], ILE)
        , ([DA, DT, DC], ILE)
        , ([DA, DT, DA], ILE)
        , ([DA, DT, DG], MET)
        , ([DA, DC, DT], THR)
        , ([DA, DC, DC], THR)
        , ([DA, DC, DA], THR)
        , ([DA, DC, DG], THR)
        , ([DA, DA, DT], ASN)
        , ([DA, DA, DC], ASN)
        , ([DA, DA, DA], LYS)
        , ([DA, DA, DG], LYS)
        , ([DA, DG, DT], SER)
        , ([DA, DG, DC], SER)
        , ([DA, DG, DA], ARG)
        , ([DA, DG, DG], ARG)
        , ([DG, DT, DT], VAL)
        , ([DG, DT, DC], VAL)
        , ([DG, DT, DA], VAL)
        , ([DG, DT, DG], VAL)
        , ([DG, DC, DT], ALA)
        , ([DG, DC, DC], ALA)
        , ([DG, DC, DA], ALA)
        , ([DG, DC, DG], ALA)
        , ([DG, DA, DT], ASP)
        , ([DG, DA, DC], ASP)
        , ([DG, DA, DA], GLU)
        , ([DG, DA, DG], GLU)
        , ([DG, DG, DT], GLY)
        , ([DG, DG, DC], GLY)
        , ([DG, DG, DA], GLY)
        , ([DG, DG, DG], GLY)
        ]

-- | taken from https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables/blob/master/codon_usage_data/tables/h_sapiens_9606.csv
codonFrequencies :: Map [DNA] Double
codonFrequencies =
    fromList
        [ ([DT, DA, DA], 0.30)
        , ([DT, DA, DG], 0.24)
        , ([DT, DG, DA], 0.47)
        , ([DG, DC, DA], 0.23)
        , ([DG, DC, DC], 0.40)
        , ([DG, DC, DG], 0.11)
        , ([DG, DC, DT], 0.27)
        , ([DT, DG, DC], 0.54)
        , ([DT, DG, DT], 0.46)
        , ([DG, DA, DC], 0.54)
        , ([DG, DA, DT], 0.46)
        , ([DG, DA, DA], 0.42)
        , ([DG, DA, DG], 0.58)
        , ([DT, DT, DC], 0.54)
        , ([DT, DT, DT], 0.46)
        , ([DG, DG, DA], 0.25)
        , ([DG, DG, DC], 0.34)
        , ([DG, DG, DG], 0.25)
        , ([DG, DG, DT], 0.16)
        , ([DC, DA, DC], 0.58)
        , ([DC, DA, DT], 0.42)
        , ([DA, DT, DA], 0.17)
        , ([DA, DT, DC], 0.47)
        , ([DA, DT, DT], 0.36)
        , ([DA, DA, DA], 0.43)
        , ([DA, DA, DG], 0.57)
        , ([DC, DT, DA], 0.07)
        , ([DC, DT, DC], 0.20)
        , ([DC, DT, DG], 0.40)
        , ([DC, DT, DT], 0.13)
        , ([DT, DT, DA], 0.08)
        , ([DT, DT, DG], 0.13)
        , ([DA, DT, DG], 1.00)
        , ([DA, DA, DC], 0.53)
        , ([DA, DA, DT], 0.47)
        , ([DC, DC, DA], 0.28)
        , ([DC, DC, DC], 0.32)
        , ([DC, DC, DG], 0.11)
        , ([DC, DC, DT], 0.29)
        , ([DC, DA, DA], 0.27)
        , ([DC, DA, DG], 0.73)
        , ([DA, DG, DA], 0.21)
        , ([DA, DG, DG], 0.21)
        , ([DC, DG, DA], 0.11)
        , ([DC, DG, DC], 0.18)
        , ([DC, DG, DG], 0.20)
        , ([DC, DG, DT], 0.08)
        , ([DA, DG, DC], 0.24)
        , ([DA, DG, DT], 0.15)
        , ([DT, DC, DA], 0.15)
        , ([DT, DC, DC], 0.22)
        , ([DT, DC, DG], 0.05)
        , ([DT, DC, DT], 0.19)
        , ([DA, DC, DA], 0.28)
        , ([DA, DC, DC], 0.36)
        , ([DA, DC, DG], 0.11)
        , ([DA, DC, DT], 0.25)
        , ([DG, DT, DA], 0.12)
        , ([DG, DT, DC], 0.24)
        , ([DG, DT, DG], 0.46)
        , ([DG, DT, DT], 0.18)
        , ([DT, DG, DG], 1.00)
        , ([DT, DA, DC], 0.56)
        , ([DT, DA, DT], 0.44)
        ]
