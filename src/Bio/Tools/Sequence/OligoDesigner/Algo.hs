module Bio.Tools.Sequence.OligoDesigner.Algo
    (split
    ) where

import           Bio.Tools.Sequence.OligoDesigner.Types (GapSize, MaxOligSize,
                                                         MinOverlap, OligSize,
                                                         OligSplitting (..),
                                                         OligsCount, Quality,
                                                         SequenceLen)
import           Debug.Trace                            (trace)


split :: SequenceLen -> MaxOligSize -> Quality -> MinOverlap -> Maybe OligSplitting
split sequLen maxOligSize quality minOverlap = do

    let minOligSize = minOverlap * 2
    let realMaxOligSize = minimum [maxOligSize, sequLen]
    let oligSizes = reverse [minOligSize .. realMaxOligSize]

    (oligSize, gapSize, count) <- findOptimalSplitting oligSizes

    let strand5 = strand5Coords oligSize gapSize [0 .. count - 1]
    let offset = oligSize - div (oligSize - gapSize) 2
    let shiftedStrand5 = shiftToOffset offset strand5
    let strand3 = init shiftedStrand5 ++ [fmap (const sequLen) (last shiftedStrand5)]

    return $ OligSplitting strand5 strand3
  where
    findOptimalSplitting :: [OligSize] -> Maybe (OligSize, GapSize, OligsCount)
    findOptimalSplitting []       = Nothing
    findOptimalSplitting (x : xs) =
        case splittingForOligSize x of
            Just (gap, count) -> return (x, gap, count)
            _                 -> findOptimalSplitting xs

    splittingForOligSize :: OligSize -> Maybe (GapSize, OligsCount)
    splittingForOligSize size = do
        let gapSizeMax = truncate $ (1 - quality) * realToFrac (size - 2 * minOverlap)
        splittingForGapSize size [0 .. gapSizeMax]

    splittingForGapSize :: OligSize -> [GapSize] -> Maybe (GapSize, OligsCount)
    splittingForGapSize _ []          = Nothing
    splittingForGapSize oligSize (x : xs) = do
            let minOligsCount = div (sequLen + x - (oligSize - minOverlap)) (oligSize + x)
            let maxOligsCount = div (sequLen + x) (oligSize + x)
            case splittingForOligsCount x oligSize [minOligsCount .. maxOligsCount] of
                Nothing           -> splittingForGapSize oligSize xs
                (Just oligsCount) -> return (x, oligsCount)

    splittingForOligsCount :: GapSize -> OligSize -> [OligsCount] -> Maybe OligsCount
    splittingForOligsCount _ _ []                               = trace "splittingForOligsCount : Nothing" Nothing
    splittingForOligsCount gapSize oligSize (x : xs) = do
        let expectedRest = oligSize - div (oligSize - gapSize) 2
        let realRest = sequLen - x * (oligSize + gapSize) + gapSize
        let leftBound = -maximum [div oligSize 10, 1]
        let rightBound = maximum [div oligSize 10, 1]
        if (realRest-expectedRest) >= leftBound && (realRest-expectedRest) <= rightBound && oligSize >= (2 * minOverlap + gapSize) && x > 0
        then return x
        else splittingForOligsCount gapSize oligSize xs


strand5Coords :: OligSize -> GapSize -> [Int] -> [(Int, Int)]
strand5Coords size gap = map toCoords where
    toCoords :: Int -> (Int, Int)
    toCoords num = (x, y) where
        x = num * (size + gap)
        y = x + size

shiftToOffset :: Int -> [(Int, Int)] -> [(Int, Int)]
shiftToOffset offset = map (\(x, y) -> (x + offset, y + offset))







