module Bio.Tools.Sequence.OligoDesigner.Algo
    (split
    ) where

import Bio.Tools.Sequence.OligoDesigner.Types (SequenceLen
                                                   ,MaxOligSize
                                                   ,MinOverlap
                                                   ,Quality
                                                   ,GapSize
                                                   ,OligsCount
                                                   ,OligSize, OligSplitting(..))
--import Debug.Trace (trace)


split :: SequenceLen -> MaxOligSize -> Quality -> MinOverlap -> Maybe OligSplitting
split sequLen maxOligSize quality minOverlap = do

    let minOligSize = minOverlap * 2
    let realMaxOligSize = minimum [maxOligSize, sequLen]
    let sizes = reverse [minOligSize .. realMaxOligSize]

    (oligSize, gapSize, count) <- findOptimalSplitting sequLen quality minOverlap sizes

    let strand5 = strand5Coords oligSize gapSize [0 .. count - 1]
    let offset = oligSize - div (oligSize - gapSize) 2
    let shiftedStrand5 = shiftToOffset offset strand5
    let strand3 = init shiftedStrand5 ++ [fmap (\_ -> sequLen) (last shiftedStrand5)]

    return $ OligSplitting strand5 strand3

strand5Coords :: OligSize -> GapSize -> [Int] -> [(Int, Int)]
strand5Coords size gap nums = map toCoords nums where
    toCoords :: Int -> (Int, Int)
    toCoords num = (x, y) where
        x = num * (size + gap)
        y = x + size

shiftToOffset :: Int -> [(Int, Int)] -> [(Int, Int)]
shiftToOffset offset coords = map (\(x, y) -> (x + offset, y + offset)) coords

findOptimalSplitting :: SequenceLen -> Quality -> MinOverlap -> [Int] -> Maybe (OligSize, GapSize, OligsCount)
findOptimalSplitting _ _ _ []                            =  Nothing
findOptimalSplitting sequLen quality minOverlap (x : xs) =  case processOligSize sequLen quality minOverlap x of
    Just (gap, count) -> return (x, gap, count)
    _                 -> findOptimalSplitting sequLen quality minOverlap xs

processOligSize :: SequenceLen -> Quality -> MinOverlap -> OligSize -> Maybe (GapSize, OligsCount)
processOligSize sequLen quality minOverlap size = do
    let gapSize = truncate $ (1 - quality) * realToFrac (size - 2 * minOverlap)
    let minOligsCount = div (sequLen + gapSize - (size - minOverlap)) (size + gapSize)
    let maxOligsCount = div (sequLen + gapSize) (size + gapSize)
    oligsCount <- countOligsCount sequLen minOverlap gapSize size [minOligsCount .. maxOligsCount]
    return (gapSize, oligsCount)

countOligsCount :: SequenceLen -> MinOverlap -> GapSize -> OligSize -> [Int] -> Maybe OligsCount
countOligsCount _ _ _ _ []                               = Nothing
countOligsCount sequLen minOverlap gapSize size (x : xs) = do
    let expectedRest = size - div (size - gapSize) 2
    let realRest = sequLen - x * (size + gapSize) + gapSize
    let leftBound = -maximum [div size 10, 1]
    let rightBound = maximum [div size 10, 1]
    if (realRest-expectedRest) >= leftBound && (realRest-expectedRest) <= rightBound && size >= (2 * minOverlap + gapSize) && x > 0
    then return x
    else countOligsCount sequLen minOverlap gapSize size xs