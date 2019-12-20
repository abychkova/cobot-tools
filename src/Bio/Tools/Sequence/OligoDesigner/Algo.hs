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
import Debug.Trace (trace)


split :: SequenceLen -> MaxOligSize -> Quality -> MinOverlap -> Maybe OligSplitting
split sequLen maxOligSize quality minOverlap = do

    let minOligSize = minOverlap * 2
    let realMaxOligSize = minimum [maxOligSize, sequLen]
    let sizes = reverse [minOligSize .. realMaxOligSize]

    (oligSize, gapSize, count) <- findOptimalSplitting sizes

    let strand5 = strand5Coords oligSize gapSize [0 .. count - 1]
    let offset = oligSize - div (oligSize - gapSize) 2
    let shiftedStrand5 = shiftToOffset offset strand5
    let strand3 = init shiftedStrand5 ++ [fmap (\_ -> sequLen) (last shiftedStrand5)]

    return $ OligSplitting strand5 strand3
  where
    findOptimalSplitting :: [Int] -> Maybe (OligSize, GapSize, OligsCount)
    findOptimalSplitting []                            =  trace ("findOptimalSplitting : Nothing") $ Nothing
    findOptimalSplitting (x : xs) =
        case processOligSize x of
            Just (gap, count) -> return (x, gap, count)
            _                 -> findOptimalSplitting xs

    processOligSize :: OligSize -> Maybe (GapSize, OligsCount)
    processOligSize size = do
        let gapSizeMax = truncate $ (1 - quality) * realToFrac (size - 2 * minOverlap)
        processGapSize size [0 .. gapSizeMax]

    processGapSize :: OligSize -> [GapSize] -> Maybe (GapSize, OligsCount)
    processGapSize _ []          = Nothing
    processGapSize size (x : xs) = do
            let minOligsCount = div (sequLen + x - (size - minOverlap)) (size + x)
            let maxOligsCount = div (sequLen + x) (size + x)
            case countOligsCount x size [minOligsCount .. maxOligsCount] of
                Nothing           -> processGapSize size xs
                (Just oligsCount) -> return (x, oligsCount)

    countOligsCount :: GapSize -> OligSize -> [Int] -> Maybe OligsCount
    countOligsCount _ _ []                               = trace ("countOligsCount : Nothing") $ Nothing
    countOligsCount gapSize size (x : xs) = do
        let expectedRest = size - div (size - gapSize) 2
        let realRest = sequLen - x * (size + gapSize) + gapSize
        let leftBound = -maximum [div size 10, 1]
        let rightBound = maximum [div size 10, 1]
        if (realRest-expectedRest) >= leftBound && (realRest-expectedRest) <= rightBound && size >= (2 * minOverlap + gapSize) && x > 0
        then return x
        else countOligsCount gapSize size xs


strand5Coords :: OligSize -> GapSize -> [Int] -> [(Int, Int)]
strand5Coords size gap nums = map toCoords nums where
    toCoords :: Int -> (Int, Int)
    toCoords num = (x, y) where
        x = num * (size + gap)
        y = x + size

shiftToOffset :: Int -> [(Int, Int)] -> [(Int, Int)]
shiftToOffset offset coords = map (\(x, y) -> (x + offset, y + offset)) coords







