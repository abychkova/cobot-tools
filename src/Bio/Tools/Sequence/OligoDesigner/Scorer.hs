module Bio.Tools.Sequence.OligoDesigner.Scorer
 (
    score
 ) where

import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligSet(..), standardTemperature)
import Bio.Tools.Sequence.ViennaRNA.Internal.Cofold (cofold)
import Debug.Trace (trace)

score :: OligSet -> Float
score (OligSet forward reversed) = score where
    all = mix forward reversed
    len = length all - 1

    diagonal = [(sequ $ all !! x, sequ $ all !! y) | x <- [0 .. len], y <- [x .. len],  x == y]
    underDiagonal = [(sequ $ all !! x, sequ $ all !! y) | x <- [0 .. len], y <- [x .. len], abs (x - y) == 1]
    matrix = [(sequ $ all !! x, sequ $ all !! y) | x <- [0 .. len], y <- [x .. len], abs (x - y) > 1]

    diagonalScore = maximum $ map (abs . fst . cofold standardTemperature) diagonal
    underDiagonalScore = minimum $ map (abs . fst . cofold standardTemperature) underDiagonal
    matrixScore = maximum $ map (abs . fst . cofold standardTemperature) matrix

    score = underDiagonalScore - max diagonalScore matrixScore
    
mix :: [a] -> [a] -> [a]
mix (x:xs) (y:ys) = x : y : mix xs ys
mix x [] = x
mix [] y = y