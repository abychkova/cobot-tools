module Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils
 (assemble
 ,buildOligSet
 ,slice
 ,translate
 ,getAAIndex
 ,mixOligs
 ,compareBySecond
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
    construct (Olig seq1 startLeft endLeft : xs) (Olig seq2 startRight endRight : ys) prevEnd acc = construct xs ys endRight res
      where
        partFormOlig1 = drop (prevEnd - startLeft) seq1
        partFormOlig2 = map cNA (drop (endLeft - startRight) (reverse seq2))
        res = acc ++ partFormOlig1 ++ partFormOlig2

buildOligSet :: OligSplitting -> [DNA] -> OligSet
buildOligSet splitting sequ = OligSet strand5' strand3' splitting
  where
    strand5' = map (buildOlig id sequ) (strand5 splitting)
    strand3' = map (buildOlig reverse (map cNA sequ)) (strand3 splitting)

    buildOlig :: ([DNA] -> [DNA]) -> [DNA] -> OligBounds -> Olig
    buildOlig fun dna (start, end) = Olig sliceDNA start end
      where
        sliceDNA = fun $ slice start end dna

--excluding end
slice :: Int -> Int -> [a] -> [a]
slice start end xs | start < 0 || end < 0 || start > end = error "incorrect coordinates"
                   | otherwise = take (end - start) (drop start xs)

translate :: [DNA] -> [DNA]
translate = map cNA
            
--TODO: test me      
getAAIndex :: Int -> Int
getAAIndex coordinate = ceiling (realToFrac (coordinate + 1) / 3)

compareBySecond :: Ord b => (a, b) -> (a, b) -> Ordering
compareBySecond p1 p2 = compare (snd p1) (snd p2)