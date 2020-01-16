{-# LANGUAGE OverloadedStrings #-}

module OligoDesigner.SpecOligoDesignerSplitter where

import           Bio.Tools.Sequence.OligoDesigner.Splitter (split)
import           Bio.Tools.Sequence.OligoDesigner.Types    (OligSplitting (..), OligSplittingConfig (..),
                                                            SequenceLen, pretty)
import           Data.Maybe                                (fromMaybe)
import           Debug.Trace                               (trace)
import           Test.Hspec                                (Spec, describe, it,
                                                            shouldBe)
import           Test.Hspec.QuickCheck                     (modifyMaxSize, prop)
import           Test.QuickCheck                           (Gen, Property,
                                                            elements, forAll)

emptySplitting :: OligSplitting
emptySplitting = OligSplitting [] []

oligoDesignerSplitterSpec :: Spec
oligoDesignerSplitterSpec =
    describe "Oligo-Designer splitter spec" $ do
        splitSequence
        splitSequence2
        splitSequence3
        splitSequence4
        splitSequence5
        realExample
        realExampleWithGap
        realExampleWithGap2
        realExampleWithGap3
        quickCheckSpec

splitSequence :: Spec
splitSequence =
    describe "splitSequence" $
    it "should correct split sequence" $ do
        let conf = OligSplittingConfig 4 1 1
        let res = fromMaybe emptySplitting (split conf 13)
        trace (pretty res) $ strand5 res `shouldBe` [(0,4), (4,8), (8,12)]
        strand3 res `shouldBe` [(2,6), (6,10), (10,13)]

splitSequence2 :: Spec
splitSequence2 =
    describe "splitSequence2" $
    it "should correct split sequence" $ do
        let conf = OligSplittingConfig 7 1 1
        let res = fromMaybe emptySplitting (split conf 13)
        trace (pretty res) $ strand5 res `shouldBe` [(0,5), (5,10)]
        strand3 res `shouldBe` [(3,8), (8,13)]

splitSequence3 :: Spec
splitSequence3 =
    describe "splitSequence3" $
    it "should correct split sequence" $ do
        let conf = OligSplittingConfig 7 1 3
        let res = split conf 13
        res `shouldBe` Nothing

splitSequence4 :: Spec
splitSequence4 =
    describe "splitSequence4" $
    it "should correct split sequence" $ do
        let conf = OligSplittingConfig 6 1 1
        let res = fromMaybe emptySplitting (split conf 11)
        trace (pretty res) $ strand5 res `shouldBe` [(0,4), (4,8)]
        strand3 res `shouldBe` [(2,6), (6,11)]

splitSequence5 :: Spec
splitSequence5 =
    describe "splitSequence5" $
    it "should correct split sequence" $ do
        let conf = OligSplittingConfig 5 1 1
        let res = fromMaybe emptySplitting (split conf 10)
        trace (pretty res) $ strand5 res `shouldBe` [(0,4), (4,8)]
        strand3 res `shouldBe` [(2,6), (6,10)]

realExample :: Spec
realExample =
    describe "realExample" $
    it "should correct split sequence" $ do
        let conf = OligSplittingConfig 60 1 18
        let res = fromMaybe emptySplitting (split conf 600)
        trace (pretty res) $ strand5 res `shouldBe` [(0,57), (57,114), (114,171), (171,228), (228,285), (285,342), (342,399), (399,456), (456,513), (513,570)]
        strand3 res `shouldBe` [(29,86), (86,143), (143,200), (200,257), (257,314), (314,371), (371,428), (428,485), (485,542), (542,600)]

realExampleWithGap :: Spec
realExampleWithGap =
    describe "realExampleWithGap" $
    it "should correct split sequence" $ do
        let conf = OligSplittingConfig 60 0.7 18
        let res = fromMaybe emptySplitting (split conf 600)
        trace (pretty res) $ strand5 res `shouldBe` [(0,60),(63,123),(126,186),(189,249),(252,312),(315,375),(378,438),(441,501),(504,564)]
        strand3 res `shouldBe` [(32,92),(95,155),(158,218),(221,281),(284,344),(347,407),(410,470),(473,533),(536,600)]

realExampleWithGap2 :: Spec
realExampleWithGap2 =
    describe "realExampleWithGap2" $
    it "should correct split sequence" $ do
        let conf = OligSplittingConfig 60 0.3 18
        let res = fromMaybe emptySplitting (split conf 600)
        trace (pretty res) $ strand5 res `shouldBe` [(0,60),(63,123),(126,186),(189,249),(252,312),(315,375),(378,438),(441,501),(504,564)]
        strand3 res `shouldBe` [(32,92),(95,155),(158,218),(221,281),(284,344),(347,407),(410,470),(473,533),(536,600)]

realExampleWithGap3 :: Spec
realExampleWithGap3 =
    describe "realExampleWithGap3" $
    it "should correct split sequence" $ do
        let conf = OligSplittingConfig 60 0.7 18
        let res = fromMaybe emptySplitting (split conf 355)
        trace (pretty res) $ strand5 res `shouldBe` [(0,60),(65,125),(130,190),(195,255),(260,320)]
        strand3 res `shouldBe` [(33,93),(98,158),(163,223),(228,288),(293,355)]

quickCheckSpec :: Spec
quickCheckSpec = modifyMaxSize (const 1000) $ prop "quickCheckSpec" splitWithRandom

splitWithRandom :: Property
splitWithRandom =
  forAll genSequenceLens
    (\seqLen -> forAll genOligMaxLens
        (\oligSize -> forAll genQualities
            (\quality -> forAll genOverlapses
                (\overlaps -> isGood overlaps seqLen (split (OligSplittingConfig oligSize (realToFrac quality / 10) overlaps) seqLen))
            )
        )
    )

isGood :: Int -> SequenceLen -> Maybe OligSplitting -> Bool
isGood _ _ Nothing  = trace "Nothing was found for this parameters" False
isGood minOverlaps sequLen (Just splitting@(OligSplitting strand5' strand3')) =
    if not res then trace (pretty splitting) res else res where
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
