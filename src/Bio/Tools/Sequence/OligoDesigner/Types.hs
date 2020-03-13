module Bio.Tools.Sequence.OligoDesigner.Types
    (SequenceLen
    ,GapSize
    ,Codon
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
    ,standardTemperature
    ,emptyMatrixCell) where

import           Bio.NucleicAcid.Nucleotide.Type      (DNA (..))
import           Bio.Tools.Sequence.CodonOptimization (CodonOptimizationConfig (..))
import           Data.List                            (foldl')
import           GHC.Generics                         (Generic)
import Control.DeepSeq (NFData)
import           Data.Default (Default (..))
import Data.Matrix (Matrix)
import Data.Text (Text)

type SequenceLen =  Int
type GapSize = Int

type OligsCount = Int
type OligSize = Int
type OligBounds = (Int, Int)
type Codon = [DNA]


--data MatrixCell = MatrixCell {olig1 :: Olig, olig2 :: Olig, rna :: Float} deriving (Show, Eq, Generic)
data MatrixCell = MatrixCell {olig1 :: OligLight, olig2 :: OligLight, rna :: Float} deriving (Show, Eq, Generic)
data OligSplitting = OligSplitting {strand5 :: [OligBounds], strand3 :: [OligBounds]} deriving (Show, Eq, NFData, Generic)
data Olig = Olig {sequDNA :: [DNA], start :: Int, end :: Int} deriving (Show, Eq, NFData, Generic)
--data OligLight = OligLight {sequStr :: Text, olig :: Olig} deriving (Show, Eq, NFData, Generic)
data OligLight = OligLight {sequStr :: String, olig :: Olig} deriving (Show, Eq, NFData, Generic)
data OligSet = OligSet {forward :: [Olig], reversed :: [Olig], coordinates :: OligSplitting} deriving (Show, Eq, NFData, Generic)

standardTemperature :: Double
standardTemperature = 37

emptyMatrixCell :: MatrixCell
--emptyMatrixCell = MatrixCell (Olig ""0 0) (Olig "" 0 0) 0
emptyMatrixCell = MatrixCell (OligLight "" (Olig ""0 0)) (OligLight "" (Olig "" 0 0)) 0

data OligsSplittingConfig = OligsSplittingConfig {
    maxOligSize    :: Int,
    overlapQuality :: Double,
    minOverlap     :: Int
}

instance Default OligsSplittingConfig where
  def = OligsSplittingConfig 60 1 18

data OligsDesignerConfig = OligsDesignerConfig {
    codonOptimizationConfig    :: CodonOptimizationConfig,
    oligSplittingConfig        :: OligsSplittingConfig,
    rnaScoreFactor             :: Double,
    oligsGCContentFactor       :: Double,
    gcContentScoreFactor       :: Double,
    maxOptimizationIteration   :: Int,
    maxFixForbiddenIteration   :: Int
}
