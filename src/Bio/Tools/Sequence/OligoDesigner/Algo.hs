module Bio.Tools.Sequence.OligoDesigner.Algo
 (
    generateOligs
 ) where

import           Bio.NucleicAcid.Nucleotide.Type                (DNA (..))
import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligSet(..), OligSplitting(..), OligBound)
import Bio.Tools.Sequence.OligoDesigner.Splitter (split)

generateOligs :: [DNA] -> Maybe OligSet
generateOligs sequ = do
    splitting <- split (length sequ) 70 1 18
    let strand5' = map (toOlig sequ) (strand5 splitting)
    let strand3' = map (toOlig sequ) (strand3 splitting)
    return $ OligSet strand5' strand3'
    
toOlig :: [DNA] -> OligBound -> Olig
toOlig sequ (start, end) = Olig oligSequ start end where
    oligSequ = slice start end sequ

slice :: Int -> Int -> [a] -> [a]
slice start end xs = take (end - start) (drop start xs)