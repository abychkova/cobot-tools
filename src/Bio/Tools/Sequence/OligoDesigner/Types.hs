module Bio.Tools.Sequence.OligoDesigner.Types
    (SequenceLen
    ,GapSize
    ,Codon
    ,MatrixCell(..)
    ,Olig(..)
    ,OligSet(..)
    ,OligsCount
    ,OligSize
    ,OligBounds
    ,OligSplitting(..)
    ,OligsDesignerConfig(..)
    ,OligsSplittingConfig(..)
    ,standardTemperature) where

import           Bio.NucleicAcid.Nucleotide.Type      (DNA (..))
import           Bio.Tools.Sequence.CodonOptimization (CodonOptimizationConfig (..))
import           Data.List                            (foldl')
import           GHC.Generics                         (Generic)
import Control.DeepSeq (NFData)
import           Data.Default (Default (..))
import Data.Matrix (Matrix)

type SequenceLen =  Int
type GapSize = Int

type OligsCount = Int
type OligSize = Int
type OligBounds = (Int, Int)
type Codon = [DNA]

data MatrixCell = MatrixCell {olig1 :: Olig, olig2 :: Olig, rna :: Float} deriving (Show, Eq, Generic)

data OligSplitting = OligSplitting {strand5 :: [OligBounds], strand3 :: [OligBounds]} deriving (Show, Eq, NFData, Generic)
data Olig = Olig {sequ :: [DNA], start :: Int, end :: Int} deriving (Show, Eq, NFData, Generic)
data OligSet = OligSet {forward :: [Olig], reversed :: [Olig], coordinates :: OligSplitting} deriving (Show, Eq, NFData, Generic)

standardTemperature :: Double
standardTemperature = 37

data OligsSplittingConfig = OligsSplittingConfig {
    maxOligSize    :: Int,
    overlapQuality :: Double,
    minOverlap     :: Int
}

instance Default OligsSplittingConfig where
  def = OligsSplittingConfig 60 1 18

data OligsDesignerConfig = OligsDesignerConfig {
    codonOptimizationConfig   :: CodonOptimizationConfig,
    oligSplittingConfig        :: OligsSplittingConfig,
    rnaScoreFactor             :: Double,
    oligsGCContentFactor       :: Double,
    gcContentScoreFactor       :: Double,
    maxOptimizationIteration   :: Int
}
