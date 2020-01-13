{-# LANGUAGE OverloadedStrings #-}

module OligoDesigner.SpecUtils where

import Bio.NucleicAcid.Nucleotide (DNA)
import Bio.Tools.Sequence.OligoDesigner.Algo (generateOligs)
import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligSet(..))
import Bio.Tools.Sequence.OligoDesigner.Utils (assemble, weightedRandom)
import Debug.Trace (trace)
import Test.Hspec (Spec, describe, it, shouldBe, shouldSatisfy)

utilsSpec :: Spec
utilsSpec =
    describe "utilsSpec" $ do
        assembleSpec
        assembleWithGapSpec
        weightedRandomSimpleSpec
        weightedRandomSpec
        weightedRandomAllZeroSpec
        weightedRandomForUnsortedSpec

assembleSpec :: Spec
assembleSpec =
    describe "assembleSpec" $
    it "assembleSpec" $ do
        let fwd = [Olig "AA" 0 2, Olig "TT" 2 4, Olig "GG" 4 6, Olig "CC" 6 8, Olig "AA" 8 10]
        let rsd = [Olig "TA" 1 3, Olig "AC" 3 5, Olig "CG" 5 7, Olig "GT" 7 9, Olig "TA" 9 11]
        let oligs = OligSet fwd rsd
        assemble oligs `shouldBe` "AATTGGCCAAT"

assembleWithGapSpec :: Spec
assembleWithGapSpec =
    describe "assembleWithGapSpec" $
    it "assembleWithGapSpec" $ do
        let fwd = [Olig "AATT" 0 4, Olig "CCAA" 6 10]
        let rsd = [Olig "ACCG" 3 7, Olig "TTA" 8 11]
        let oligs = OligSet fwd rsd
        assemble oligs `shouldBe` "AATTGGCCAAT"

weightedRandomSimpleSpec :: Spec
weightedRandomSimpleSpec =
    describe "weightedRandomSimpleSpec" $
    it "should return weigthed random element" $ do
        res <- weightedRandom [(1, 1), (2, 0)]
        res `shouldBe` 1
        res2 <- weightedRandom [(1, 0.99999999), (2, 0.00000001)]
        res2 `shouldBe` 1
        res3 <- weightedRandom [(1, 0), (2, 1)]
        res3 `shouldBe` 2

weightedRandomSpec :: Spec
weightedRandomSpec =
    describe "weightedRandomSpec" $
    it "should return weigthed random element" $ do
        res <- sequence [weightedRandom [(3, 7), (4, 300), (1, 1600)] | x <- [0..10]]
        let count1 = length $ filter (== 1) res
        let count3 = length $ filter (== 3) res
        let count4 = length $ filter (== 4) res
        count1 `shouldSatisfy` (> count3)
        count1 `shouldSatisfy` (> count4)
        count4 `shouldSatisfy` (>= count3)

weightedRandomForUnsortedSpec :: Spec
weightedRandomForUnsortedSpec =
    describe "weightedRandomForUnsortedSpec" $
    it "should return weigthed random element" $ do
        res <- sequence [weightedRandom [(4, 300), (3, 7), (1, 1600)] | x <- [0..10]]
        let count1 = length $ filter (== 1) res
        let count3 = length $ filter (== 3) res
        let count4 = length $ filter (== 4) res
        count1 `shouldSatisfy` (> count3)
        count1 `shouldSatisfy` (> count4)
        count4 `shouldSatisfy` (>= count3)

weightedRandomAllZeroSpec :: Spec
weightedRandomAllZeroSpec =
    describe "weightedRandomAllZeroSpec" $
    it "should return weigthed random element" $ do
        res <- weightedRandom [(1, 0), (2, 0)]
        res `shouldBe` 1