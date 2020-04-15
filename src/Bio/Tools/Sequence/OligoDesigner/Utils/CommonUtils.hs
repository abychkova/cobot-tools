module Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils
    ( assemble
    , buildOligSet
    , slice
    , translate
    , getAAIndex
    , mixOligs
    , compareBySecond
    , orderByScore
    , isEqual
    ) where

import Control.Monad.Except (Except, throwError)
import Data.List            (maximumBy, minimumBy)
import Data.Ord             (comparing)

import Bio.NucleicAcid.Nucleotide.Type        (DNA (..), cNA)
import Bio.Tools.Sequence.OligoDesigner.Types (Olig (..), OligBounds, OligSet (..),
                                               OligSplitting (..), OligoDesignerError(..))


mixOligs :: OligSet -> [Olig]
mixOligs (OligSet forward reversed _) = mix forward reversed
  where
    mix :: [a] -> [a] -> [a]
    mix (x:xs) (y:ys) = x : y : mix xs ys
    mix x []          = x
    mix [] y          = y

assemble :: OligSet -> [DNA]
assemble (OligSet fwd rvd _) = construct fwd rvd 0 [] where

    construct :: [Olig] -> [Olig] -> Int -> [DNA] -> [DNA]
    construct _ [] _ acc = acc
    construct [] _ _ acc = acc
    construct (oligLeft : xs) (oligRight : ys) prevEnd acc = construct xs ys endRight res
      where
        (Olig seqLeft startLeft endLeft) = oligLeft
        (Olig seqRight startRight endRight) = oligRight
        partFormOlig1 = drop (prevEnd - startLeft) seqLeft
        partFormOlig2 = map cNA (drop (endLeft - startRight) (reverse seqRight))
        res = acc ++ partFormOlig1 ++ partFormOlig2

buildOligSet :: OligSplitting -> [DNA] -> Except OligoDesignerError OligSet
buildOligSet splitting sequ = do
    strand5' <- mapM (buildOlig id sequ) (strand5 splitting)
    strand3' <- mapM (buildOlig reverse (map cNA sequ)) (strand3 splitting)
    return $ OligSet strand5' strand3' splitting
  where
    buildOlig :: ([DNA] -> [DNA]) -> [DNA] -> OligBounds -> Except OligoDesignerError Olig
    buildOlig fun dna (start, end) = do
        sliceDNA <- slice start end dna
        return $ Olig (fun sliceDNA) start end

--excluding end
slice :: Int -> Int -> [a] -> Except OligoDesignerError [a]
slice start end xs | start < 0 || end < 0 || start > end = throwError (InvalidInterval (start, end))
                   | otherwise = return $ take (end - start) (drop start xs)

translate :: [DNA] -> [DNA]
translate = map cNA

getAAIndex :: Int -> Except OligoDesignerError Int
getAAIndex coordinate | coordinate < 0 = throwError (InvalidInterval (coordinate, coordinate))
                      | otherwise      = return $ ceiling (realToFrac (coordinate + 1) / 3 :: Double)


compareBySecond :: Ord b => (a, b) -> (a, b) -> Ordering
compareBySecond = comparing snd

orderByScore :: Ord b => [a] -> (a -> b) -> (a, a)
orderByScore sequ func = (fst $ minimumBy compareBySecond pairs, fst $ maximumBy compareBySecond pairs)
  where
    pairs = map (\value -> (value, func value)) sequ
    
isEqual :: Double -> Double -> Bool
isEqual a b = abs(a - b) < 0.00001
