module OligoDesigner.SpecSplitter (
    oligoDesignerSplitterSpec
) where

import Data.Maybe            (fromMaybe)
import Test.Hspec            (Spec, describe, it, shouldBe)
import Test.Hspec.QuickCheck (modifyMaxSize, prop)
import Test.QuickCheck       (Gen, Property, elements, forAll)

import Bio.Tools.Sequence.OligoDesigner.Splitter         (split)
import Bio.Tools.Sequence.OligoDesigner.Types            (OligSplitting (..),
                                                          OligsSplittingConfig (..), SequenceLen)

emptySplitting :: OligSplitting
emptySplitting = OligSplitting [] []

oligoDesignerSplitterSpec :: Spec
oligoDesignerSplitterSpec =
    describe "Oligo-Designer splitter spec" $ do
        splitSequence
        realExampleWithGap
        quickCheckSpec

splitSequence :: Spec
splitSequence =
    describe "splitSequence" $
    it "should correct split sequence for 4 1 1"
        (do let conf = OligsSplittingConfig 4 1 1
            let res = fromMaybe emptySplitting (split conf 13)
            strand5 res `shouldBe` [(0, 4), (4, 8), (8, 12)]
            strand3 res `shouldBe` [(2, 6), (6, 10), (10, 13)]) >>
    it "should correct split sequence for 7 1 1"
        (do let conf = OligsSplittingConfig 7 1 1
            let res = fromMaybe emptySplitting (split conf 13)
            strand5 res `shouldBe` [(0, 5), (5, 10)]
            strand3 res `shouldBe` [(3, 8), (8, 13)]) >>
    it "should correct split sequence for 7 1 3"
        (do let conf = OligsSplittingConfig 7 1 3
            let res = split conf 13
            res `shouldBe` Nothing) >>
    it "should correct split sequence for 6 1 1"
        (do let conf = OligsSplittingConfig 6 1 1
            let res = fromMaybe emptySplitting (split conf 11)
            strand5 res `shouldBe` [(0, 4), (4, 8)]
            strand3 res `shouldBe` [(2, 6), (6, 11)]) >>
    it "should correct split sequence for 5 1 1"
        (do let conf = OligsSplittingConfig 5 1 1
            let res = fromMaybe emptySplitting (split conf 10)
            strand5 res `shouldBe` [(0, 4), (4, 8)]
            strand3 res `shouldBe` [(2, 6), (6, 10)])

realExampleWithGap :: Spec
realExampleWithGap =
    describe "realExample" $
    it "should correct split sequence"
        (do let conf = OligsSplittingConfig 60 1 18
            let res = fromMaybe emptySplitting (split conf 600)
            strand5 res `shouldBe` [(0,57), (57,114), (114,171), (171,228), (228,285), (285,342), (342,399), (399,456), (456,513), (513,570)]
            strand3 res `shouldBe` [(29,86), (86,143), (143,200), (200,257), (257,314), (314,371), (371,428), (428,485), (485,542), (542,600)]) >>
    it "should correct split sequence with quality 0.7"
        (do let conf = OligsSplittingConfig 60 0.7 18
            let res = fromMaybe emptySplitting (split conf 600)
            strand5 res `shouldBe` [(0,60),(63,123),(126,186),(189,249),(252,312),(315,375),(378,438),(441,501),(504,564)]
            strand3 res `shouldBe` [(32,92),(95,155),(158,218),(221,281),(284,344),(347,407),(410,470),(473,533),(536,600)]) >>
    it "should correct split sequence with quality 0.3"
        (do let conf = OligsSplittingConfig 60 0.3 18
            let res = fromMaybe emptySplitting (split conf 600)
            strand5 res `shouldBe` [(0,60),(63,123),(126,186),(189,249),(252,312),(315,375),(378,438),(441,501),(504,564)]
            strand3 res `shouldBe` [(32,92),(95,155),(158,218),(221,281),(284,344),(347,407),(410,470),(473,533),(536,600)]) >>
    it "should correct split sequence with quality 0.7 and len 355"
        (do let conf = OligsSplittingConfig 60 0.7 18
            let res = fromMaybe emptySplitting (split conf 355)
            strand5 res `shouldBe` [(0,60),(65,125),(130,190),(195,255),(260,320)]
            strand3 res `shouldBe` [(33,93),(98,158),(163,223),(228,288),(293,355)])

quickCheckSpec :: Spec
quickCheckSpec = modifyMaxSize (const 1000) $ prop "quickCheckSpec" splitWithRandom

splitWithRandom :: Property
splitWithRandom =
  forAll genSequenceLens
    (\seqLen -> forAll genOligMaxLens
        (\oligSize -> forAll genQualities
            (\quality -> forAll genOverlapses
                (\overlaps -> isGood overlaps seqLen (split (OligsSplittingConfig oligSize (realToFrac quality / 10) overlaps) seqLen))
            )
        )
    )

isGood :: Int -> SequenceLen -> Maybe OligSplitting -> Bool
isGood _ _ Nothing  = False
isGood minOverlaps sequLen (Just (OligSplitting strand5' strand3')) =
    if not res then res else res where
        zip53 = zip strand5' strand3'
        zip35 = zip strand3' (tail strand5')

        overlaps :: [((Int, Int), (Int, Int))] -> [Int]
        overlaps = map (\(olig1, olig2) -> snd olig1 - fst olig2)

        overlaps53 = overlaps zip53
        overlaps35 = overlaps zip35

        zipOverlaps = zip overlaps53 overlaps35

        res = all (>= minOverlaps) overlaps53 &&
                all (>= minOverlaps) overlaps35 &&
                all (\(x, y) -> x - y < 1) zipOverlaps &&
                sequLen == snd (last strand3')

qualities :: [Int]
qualities = [0..10]

sequenceLens :: [Int]
sequenceLens = [200..2000]

oligMaxLens :: [Int]
oligMaxLens = [70..200]

overlapses :: [Int]
overlapses = [18..30]

genQualities :: Gen Int
genQualities = elements qualities

genOverlapses :: Gen Int
genOverlapses = elements overlapses

genOligMaxLens :: Gen Int
genOligMaxLens = elements oligMaxLens

genSequenceLens :: Gen Int
genSequenceLens = elements sequenceLens
