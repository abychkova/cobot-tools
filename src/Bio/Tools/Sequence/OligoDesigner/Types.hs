module Bio.Tools.Sequence.OligoDesigner.Types
    (SequenceLen
    ,MaxOligSize
    ,MinOverlap
    ,Quality
    ,GapSize
    ,Olig(..)
    ,OligSet(..)
    ,OligsCount
    ,OligSize
    ,OligBound
    ,OligSplitting(..)
    ,standardTemperature
    ,pretty) where

import           GHC.Generics                (Generic)
import Data.List (foldl')
import Debug.Trace (trace)
import           Bio.NucleicAcid.Nucleotide.Type                (DNA (..))

type SequenceLen =  Int
type MaxOligSize =  Int
type MinOverlap =  Int
type Quality =  Double
type GapSize = Int

type OligsCount = Int
type OligSize = Int
type OligBound = (Int, Int)

data OligSplitting = OligSplitting {strand5 :: [OligBound], strand3 :: [OligBound]} deriving (Show, Eq, Generic)

pretty :: OligSplitting -> String
pretty (OligSplitting strand5 strand3) = "5': " ++ str5 ++ "\n3': " ++ str3 where
    str5 = fst $ foldl' concat ("", 0) strand5
    str3 = fst $ foldl' concat ("", 0) strand3

    concat :: (String, Int) -> (Int, Int) -> (String, Int)
    concat (res, prev) p@(x, y) = (res ++ replicate (x - prev) ' ' ++ "(" ++ replicate (y - x) '_' ++ ")", y)

data Olig = Olig {sequ :: [DNA], start :: Int, end :: Int} deriving (Show, Eq, Generic)
data OligSet = OligSet {forward :: [Olig], reversed :: [Olig]}

standardTemperature :: Double
standardTemperature = 37