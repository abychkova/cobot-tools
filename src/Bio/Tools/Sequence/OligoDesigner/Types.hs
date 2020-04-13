module Bio.Tools.Sequence.OligoDesigner.Types
    (SequenceLen
    ,GapSize
    ,Codon
    ,TargetGC
    ,Weight
    ,MatrixCell(..)
    ,Olig(..)
    ,OligLight(..)
    ,OligSet(..)
    ,OligsCount
    ,OligSize
    ,OligBounds
    ,OligSplitting(..)
    ,OligsDesignerConfig(..)
    ,OligsSplittingConfig(..)
    ,OligsDesignerInnerConfig(..)
    ,standardTemperature
    ,emptyMatrixCell) where

import           Bio.NucleicAcid.Nucleotide.Type            (DNA (..))
import           Bio.Tools.Sequence.CodonOptimization       (CodonOptimizationConfig (..))
import           Bio.Tools.Sequence.CodonOptimization.Types (Organism)
import           Control.DeepSeq                            (NFData)
import           Data.Default                               (Default (..))
import           GHC.Generics                               (Generic)
import           Text.Regex.TDFA                            (Regex)

type SequenceLen = Int
type GapSize = Int

type TargetGC = Double
type Weight = Double
type OligsCount = Int
type OligSize = Int
type OligBounds = (Int, Int)
type Codon = [DNA]

data MatrixCell =
    MatrixCell
        { olig1 :: OligLight
        , olig2 :: OligLight
        , rna :: Float
        }
    deriving (Show, Eq, Generic)

data OligSplitting =
    OligSplitting
        { strand5 :: [OligBounds]
        , strand3 :: [OligBounds]
        }
    deriving (Show, Eq, NFData, Generic)

data Olig =
    Olig
        { sequDNA :: [DNA]
        , start :: Int
        , end :: Int
        }
    deriving (Show, Eq, NFData, Generic)

data OligLight =
    OligLight
        { sequStr :: String
        , olig :: Olig
        }
    deriving (Show, Eq, NFData, Generic)

instance Default OligLight where
    def = OligLight "" (Olig "" 0 0)

data OligSet =
    OligSet
        { forward :: [Olig]
        , reversed :: [Olig]
        , coordinates :: OligSplitting
        }
    deriving (Show, Eq, NFData, Generic)

standardTemperature :: Double
standardTemperature = 37

emptyMatrixCell :: MatrixCell
emptyMatrixCell = MatrixCell def def 0

data OligsSplittingConfig =
    OligsSplittingConfig
        { maxOligSize :: Int
        , overlapQuality :: Double
        , minOverlap :: Int
        }

instance Default OligsSplittingConfig where
    def = OligsSplittingConfig 60 1 18

data OligsDesignerConfig =
    OligsDesignerConfig
        { codonOptimizationConfig :: CodonOptimizationConfig
        , oligSplittingConfig :: OligsSplittingConfig
        , maxOptimizationIteration :: Int
        , maxFixForbiddenIteration :: Int
        }

data OligsDesignerInnerConfig =
    OligsDesignerInnerConfig
        { organismType :: Organism
        , gcTarget :: Double
        , regexes :: [Regex]
        , optimizationIteration :: Int
        , fixForbiddenIteration :: Int
        }