{-# LANGUAGE OverloadedStrings #-}

module SpecOligoDesigner where

import           Bio.Tools.Sequence.OligoDesigner.Types         (OligSplitting(..))
import           Bio.Tools.Sequence.OligoDesigner.Algo          (split)
import           Test.Hspec                                     (Expectation,
                                                                 Spec, describe,
                                                                 it, shouldBe,
                                                                 shouldSatisfy)
import Data.Maybe (fromMaybe, isJust)
import Debug.Trace (trace, traceShow)
import Test.QuickCheck
import Test.Hspec.QuickCheck (prop)

emptySplitting :: OligSplitting
emptySplitting = OligSplitting [] []

oligoDesignerSpec :: Spec
oligoDesignerSpec =
    describe "Oligo-Designer spec" $ do
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
        let res = fromMaybe emptySplitting (split 13 4 1 1)
        trace (show res) $ strand5 res `shouldBe` [(0,4), (4,8), (8,12)]
        strand3 res `shouldBe` [(2,6), (6,10), (10,13)]

splitSequence2 :: Spec
splitSequence2 =
    describe "splitSequence2" $
    it "should correct split sequence" $ do
        let res = fromMaybe emptySplitting (split 13 7 1 1)
        trace (show res) $ strand5 res `shouldBe` [(0,5), (5,10)]
        strand3 res `shouldBe` [(3,8), (8,13)]

splitSequence3 :: Spec
splitSequence3 =
    describe "splitSequence3" $
    it "should correct split sequence" $ do
        let res = split 13 7 1 3
        res `shouldBe` Nothing

splitSequence4 :: Spec
splitSequence4 =
    describe "splitSequence4" $
    it "should correct split sequence" $ do
        let res = fromMaybe emptySplitting (split 11 6 1 1)
        trace (show res) $ strand5 res `shouldBe` [(0,4), (4,8)]
        strand3 res `shouldBe` [(2,6), (6,11)]

splitSequence5 :: Spec
splitSequence5 =
    describe "splitSequence5" $
    it "should correct split sequence" $ do
        let res = fromMaybe emptySplitting (split 10 5 1 1)
        trace (show res) $ strand5 res `shouldBe` [(0,4), (4,8)]
        strand3 res `shouldBe` [(2,6), (6,10)]

realExample :: Spec
realExample =
    describe "realExample" $
    it "should correct split sequence" $ do
        let res = fromMaybe emptySplitting (split 600 60 1 18)
        trace (show res) $ strand5 res `shouldBe` [(0,57), (57,114), (114,171), (171,228), (228,285), (285,342), (342,399), (399,456), (456,513), (513,570)]
        strand3 res `shouldBe` [(29,86), (86,143), (143,200), (200,257), (257,314), (314,371), (371,428), (428,485), (485,542), (542,600)]

realExampleWithGap :: Spec
realExampleWithGap =
    describe "realExampleWithGap" $
    it "should correct split sequence" $ do
        let res = fromMaybe emptySplitting (split 600 60 0.7 18)
        trace (show res) $ strand5 res `shouldBe` [(0,58), (64,122), (128,186), (192,250), (256,314), (320,378), (384,442), (448,506), (512,570)]
        strand3 res `shouldBe` [(32,90), (96,154), (160,218), (224,282), (288,346), (352,410), (416,474), (480,538), (544,600)]

realExampleWithGap2 :: Spec
realExampleWithGap2 =
    describe "realExampleWithGap2" $
    it "should correct split sequence" $ do
        let res = fromMaybe emptySplitting (split 600 60 0.3 18)
        trace (show res) $ strand5 res `shouldBe` [(0,53), (64,117), (128,181), (192,245), (256,309), (320,373), (384,437), (448,501), (512,565)]
        strand3 res `shouldBe` [(32,85), (96,149), (160,213), (224,277), (288,341), (352,405), (416,469), (480,533), (544,600)]

realExampleWithGap3 :: Spec
realExampleWithGap3 =
    describe "realExampleWithGap3" $
    it "should correct split sequence" $ do
        let res = fromMaybe emptySplitting (split 355 60 0.7 18)
        trace (show res) $ strand5 res `shouldBe` [(0,59), (65,124), (130,189), (195,254), (260,319)]
        strand3 res `shouldBe` [(33,92), (98,157), (163,222), (228,287), (293,355)]

quickCheckSpec :: Spec
quickCheckSpec = prop "quickCheckSpec" $ splitWithRandom

splitWithRandom :: Property
splitWithRandom =
  forAll genSequenceLens
    (\seqLen -> forAll genOligMaxLens
        (\oligSize -> forAll genOverlapses
            (\overlaps -> forAll genQualities
                (\quality -> isJust $ split seqLen oligSize 1 overlaps)
            )
        )
    )

qualities :: [Int]
qualities = [0,1..1]

sequenceLens :: [Int]
sequenceLens = [100..1000]

oligMaxLens :: [Int]
oligMaxLens = [40..100]

overlapses :: [Int]
overlapses = [3..20]

genQualities :: Gen Int
genQualities = elements qualities

genOverlapses :: Gen Int
genOverlapses = elements overlapses

genOligMaxLens :: Gen Int
genOligMaxLens = elements oligMaxLens

genSequenceLens :: Gen Int
genSequenceLens = elements sequenceLens