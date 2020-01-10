module Bio.Tools.Sequence.OligoDesigner.Algo
 (generateOligs
 ) where

import           Bio.NucleicAcid.Nucleotide.Type                (DNA (..))
import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligSet(..), OligSplitting(..), OligBounds)
import Bio.Tools.Sequence.OligoDesigner.Utils (translateDNA, assemble)
import Bio.Tools.Sequence.OligoDesigner.Splitter (split)
import Bio.Tools.Sequence.OligoDesigner.Optimizer (optimize)
import Bio.Tools.Sequence.CodonOptimization.Algo (optimizeAA, scoreSequence)
import Bio.Protein.AminoAcid (AA)
import Bio.Tools.Sequence.CodonOptimization.Types (CodonScoreConfig)
import Bio.Tools.Sequence.OligoDesigner.Scorer (rnaCofoldScore, gcScore)
import GHC.Float (float2Double)
import Data.List (maximumBy)
import Data.Default (def)

generateOligs :: [DNA] -> Maybe OligSet
generateOligs sequ = do
    splitting <- split (length sequ) 70 1 18
    let strand5' = map (toOlig sequ) (strand5 splitting)
    let strand3' = map (toOlig $ translateDNA sequ) (strand3 splitting)
    return $ OligSet strand5' strand3'
  where
    toOlig :: [DNA] -> OligBounds -> Olig
    toOlig dna (start, end) = Olig (slice start end dna) start end

    slice :: Int -> Int -> [a] -> [a]
    slice start end xs = take (end - start) (drop start xs)

generateAndOptimizeOligs :: [AA] -> Maybe OligSet
generateAndOptimizeOligs aa = do
    let dna = optimizeAA def aa
    oligs <- generateOligs dna
    return $ maximumBy scoreCmp [optimize oligs | _ <- [0..5]] --TODO: better do this with threads
  where
    commonScore :: OligSet -> Double --TODO: better move to Scorer.hs as globalScore and add codon-optimization-score
    commonScore oligs = if rna <= 0 then realToFrac rna else gc ** 50 * rna
      where
        rna = realToFrac $ rnaCofoldScore oligs
        gc = gcScore oligs 60
        codonOpt = scoreSequence $ def assemble oligs

    scoreCmp :: OligSet -> OligSet -> Ordering
    scoreCmp oligs1 oligs2 = compare score1 score2
      where
        score1 = commonScore oligs1
        score2 = commonScore oligs2