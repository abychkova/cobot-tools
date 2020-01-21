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
    ,OligoDesignerConfig(..)
    ,OligSplittingConfig(..)
    ,standardTemperature
    ,pretty) where

import           Bio.NucleicAcid.Nucleotide.Type      (DNA (..))
import           Bio.Tools.Sequence.CodonOptimization (CodonOptimizationConfig (..))
import           Data.List                            (foldl')
import           GHC.Generics                         (Generic)
import Control.DeepSeq (NFData)

type SequenceLen =  Int
type GapSize = Int

type OligsCount = Int
type OligSize = Int
type OligBounds = (Int, Int)
type Codon = [DNA]

data MatrixCell = MatrixCell {olig1 :: Olig, olig2 :: Olig, rna :: Float} deriving (Show, Eq, Generic)

data OligSplitting = OligSplitting {strand5 :: [OligBounds], strand3 :: [OligBounds]} deriving (Show, Eq, NFData, Generic)

pretty :: OligSplitting -> String
pretty (OligSplitting strand5 strand3) = "5': " ++ str5 ++ "\n3': " ++ str3 where
    str5 = fst $ foldl' conc ("", 0) strand5
    str3 = fst $ foldl' conc ("", 0) strand3

    conc :: (String, Int) -> (Int, Int) -> (String, Int)
    conc (res, prev) (x, y) = (res ++ replicate (x - prev) ' ' ++ "(" ++ replicate (y - x) '_' ++ ")", y)

data Olig = Olig {sequ :: [DNA], start :: Int, end :: Int} deriving (Show, Eq, NFData, Generic)
data OligSet = OligSet {forward :: [Olig], reversed :: [Olig], coordinates :: OligSplitting} deriving (Show, Eq, NFData, Generic)

standardTemperature :: Double
standardTemperature = 37

data OligSplittingConfig = OligSplittingConfig {
    maxOligSize    :: Int,
    overlapQuality :: Double,
    minOverlap     :: Int
}

data OligoDesignerConfig = OligoDesignerConfig {
    codonOptimizationConfing :: CodonOptimizationConfig,
    balanceFactor            :: Double,
    oligSplittingConfig      :: OligSplittingConfig
}
