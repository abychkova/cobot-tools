module Bio.Tools.Sequence.OligoDesigner.Algo
 (generateOligs
 ) where

import           Bio.NucleicAcid.Nucleotide.Type                (DNA (..))
import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligSet(..), OligSplitting(..), OligBounds)
import Bio.Tools.Sequence.OligoDesigner.Utils (translate)
import Bio.Tools.Sequence.OligoDesigner.Splitter (split)

generateOligs :: [DNA] -> Maybe OligSet
generateOligs sequ = do
    splitting <- split (length sequ) 70 1 18
    let translated = translate sequ
    let strand5' = map (toOlig sequ) (strand5 splitting)
    let strand3' = map (toOlig translated) (strand3 splitting)
    return $ OligSet strand5' strand3'
  where
    toOlig :: [DNA] -> OligBounds -> Olig
    toOlig dna (start, end) = Olig (slice start end dna) start end

    slice :: Int -> Int -> [a] -> [a]
    slice start end xs = take (end - start) (drop start xs)