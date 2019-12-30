module Bio.Tools.Sequence.OligoDesigner.Scorer
 (rnaCofoldScore
 ,gcScore
 ) where

import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligSet(..), standardTemperature)
import Bio.Tools.Sequence.ViennaRNA.Internal.Cofold (cofold)
import Bio.NucleicAcid.Nucleotide.Type (DNA(..))

rnaCofoldScore :: OligSet -> Float
rnaCofoldScore (OligSet forward reversed) = score where
    allOligs = mix forward reversed
    len = length allOligs - 1

    diagonal      = [(sequ $ allOligs !! x, sequ $ allOligs !! y) | x <- [0 .. len], y <- [x .. len],  x == y]
    underDiagonal = [(sequ $ allOligs !! x, sequ $ allOligs !! y) | x <- [0 .. len], y <- [x .. len], abs (x - y) == 1]
    matrix        = [(sequ $ allOligs !! x, sequ $ allOligs !! y) | x <- [0 .. len], y <- [x .. len], abs (x - y) > 1]

    diagonalScore      = maximum $ map (abs . fst . cofold standardTemperature) diagonal
    underDiagonalScore = minimum $ map (abs . fst . cofold standardTemperature) underDiagonal
    matrixScore        = maximum $ map (abs . fst . cofold standardTemperature) matrix

    score = underDiagonalScore - max diagonalScore matrixScore
    
    mix :: [a] -> [a] -> [a]
    mix (x:xs) (y:ys) = x : y : mix xs ys
    mix x [] = x
    mix [] y = y

gcScore :: [DNA] -> Double -> Double
gcScore sequ target = score where
    gc = length $ filter (\nk -> nk == DC || nk == DG) sequ
    gcContent = realToFrac gc / realToFrac (length sequ)
    deltaGC = abs(target - gcContent)
    score = 1 - deltaGC / target