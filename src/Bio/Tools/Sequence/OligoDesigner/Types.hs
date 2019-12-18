module Bio.Tools.Sequence.OligoDesigner.Types
    (SequenceLen
    ,MaxOligSize
    ,MinOverlap
    ,Quality
    ,GapSize
    ,OligsCount
    ,OligSize
    ,Olig
    ,OligSplitting(..)) where

import           GHC.Generics                (Generic)
import Data.List (foldl')
import Debug.Trace (trace)

type SequenceLen =  Int
type MaxOligSize =  Int
type MinOverlap =  Int
type Quality =  Double
type GapSize = Int

type OligsCount = Int
type OligSize = Int

type Olig = (Int, Int)
data OligSplitting = OligSplitting {strand5 :: [Olig], strand3 :: [Olig]} deriving (Eq, Generic)

instance Show OligSplitting where
    show (OligSplitting strand5 strand3) = "5': " ++ str5 ++ "\n3': " ++ str3 where
        str5 = fst $ foldl' concat ("", 0) strand5
        str3 = fst $ foldl' concat ("", 0) strand3

        concat :: (String, Int) -> (Int, Int) -> (String, Int)
        concat (res, prev) p@(x, y) = (res ++ replicate (x - prev) ' ' ++ "(" ++ replicate (y - x) '_' ++ ")", y)
