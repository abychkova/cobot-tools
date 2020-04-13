module Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils
 (assemble
 ,buildOligSet
 ,slice
 ,translate
 ,getAAIndex
 ,mixOligs
 ,compareBySecond
 ,orderByScore
 ) where

import           Bio.NucleicAcid.Nucleotide.Type                (DNA (..), cNA)
import           Bio.Protein.AminoAcid.Type                     (AA (..))
import           Bio.Tools.Sequence.CodonOptimization.Constants (ak2Codon,
                                                                 codon2ak,
                                                                 codonFrequencies)
import           Bio.Tools.Sequence.CodonOptimization.Types     (Organism)
import           Bio.Tools.Sequence.OligoDesigner.Types         (Codon,
                                                                 Olig (..),
                                                                 OligBounds,
                                                                 OligSet (..),
                                                                 OligSplitting (..), OligsDesignerConfig)
import           Control.Monad.State                            (State, get,
                                                                 put)
import           Data.List                                      (sortOn, nub)
import           Data.Map                                       as Map (lookup)
import           Data.Maybe                                     (fromMaybe)
import           System.Random                                  (StdGen,
                                                                 randomR)
import Debug.Trace 
import Bio.Tools.Sequence.CodonOptimization (CodonOptimizationConfig(..))
import Text.Regex.TDFA (Regex, makeRegex, match)
import Control.Monad.State.Lazy (gets)
import System.Random.Shuffle (shuffle')
import Control.Monad.Except (Except, throwError)
import Data.List (maximumBy, minimumBy, findIndex)
        
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
    construct (Olig seqLeft startLeft endLeft : xs) (Olig seqRight startRight endRight : ys) prevEnd acc = construct xs ys endRight res
      where
        partFormOlig1 = drop (prevEnd - startLeft) seqLeft
        partFormOlig2 = map cNA (drop (endLeft - startRight) (reverse seqRight))
        res = acc ++ partFormOlig1 ++ partFormOlig2

buildOligSet :: OligSplitting -> [DNA] -> Except String OligSet
buildOligSet splitting sequ = do
    strand5' <- mapM (buildOlig id sequ) (strand5 splitting)
    strand3' <- mapM (buildOlig reverse (map cNA sequ)) (strand3 splitting)
    return $ OligSet strand5' strand3' splitting
  where
    buildOlig :: ([DNA] -> [DNA]) -> [DNA] -> OligBounds -> Except String Olig
    buildOlig fun dna (start, end) = do
        sliceDNA <- slice start end dna
        return $ Olig (fun sliceDNA) start end

--excluding end
slice :: Int -> Int -> [a] -> Except String [a]
slice start end xs | start < 0 || end < 0 || start > end = throwError "incorrect coordinates"
                   | otherwise = return $ take (end - start) (drop start xs)

translate :: [DNA] -> [DNA]
translate = map cNA
              
getAAIndex :: Int -> Except String Int
getAAIndex coordinate | coordinate < 0 = throwError "incorrect coordinates"
                      | otherwise      = return $ ceiling (realToFrac (coordinate + 1) / 3)
                      

compareBySecond :: Ord b => (a, b) -> (a, b) -> Ordering
compareBySecond p1 p2 = compare (snd p1) (snd p2)

orderByScore :: Ord b => [a] -> (a -> b) -> (a, a)
orderByScore sequence func = (fst $ minimumBy compareBySecond pairs, fst $ maximumBy compareBySecond pairs)
  where
    pairs = map (\value -> (value, func value)) sequence