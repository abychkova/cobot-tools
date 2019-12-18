{-# LANGUAGE OverloadedStrings #-}

module SpecOligoDesigner where

import           Bio.Tools.Sequence.OligoDesigner.Algo          (split, OligSplitting(..))
import           Test.Hspec                                     (Expectation,
                                                                 Spec, describe,
                                                                 it, shouldBe,
                                                                 shouldSatisfy)
import Data.Maybe (fromMaybe)

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

splitSequence :: Spec
splitSequence =
    describe "splitSequence" $
    it "should correct split sequence" $ do
        let res = fromMaybe emptySplitting (split 13 4 1 1)
        strand5 res `shouldBe` [(0,4), (4,8), (8,12)]
        strand3 res `shouldBe` [(2,6), (6,10), (10,13)]

splitSequence2 :: Spec
splitSequence2 =
    describe "splitSequence2" $
    it "should correct split sequence" $ do
        let res = fromMaybe emptySplitting (split 13 7 1 1)
        strand5 res `shouldBe` [(0,5), (5,10)]
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
        strand5 res `shouldBe` [(0,4), (4,8)]
        strand3 res `shouldBe` [(2,6), (6,11)]

splitSequence5 :: Spec
splitSequence5 =
    describe "splitSequence5" $
    it "should correct split sequence" $ do
        let res = fromMaybe emptySplitting (split 10 5 1 1)
        strand5 res `shouldBe` [(0,4), (4,8)]
        strand3 res `shouldBe` [(2,6), (6,10)]

realExample :: Spec
realExample =
    describe "realExample" $
    it "should correct split sequence" $ do
        let res = fromMaybe emptySplitting (split 600 60 1 18)
        strand5 res `shouldBe` [(0,57), (57,114), (114,171), (171,228), (228,285), (285,342), (342,399), (399,456), (456,513), (513,570)]
        strand3 res `shouldBe` [(29,86), (86,143), (143,200), (200,257), (257,314), (314,371), (371,428), (428,485), (485,542), (542,600)]