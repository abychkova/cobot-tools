{-# LANGUAGE OverloadedStrings #-}

module OligoDesigner.SpecUtils (
    utilsSpec
) where

import           Bio.Tools.Sequence.OligoDesigner.Types (Olig (..),
                                                         OligSet (..),
                                                         OligSplitting (..))
import           Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils (assemble,
                                                         slice,
                                                         buildOligSet,translate)
import           Control.Exception                      (evaluate)
import           Control.Monad.State                    (State, evalState, get,
                                                         put, runState)
import           Data.List                              (nub)
import           System.Random                          (StdGen, getStdGen, mkStdGen)
import           Test.Hspec                             (Spec, describe,
                                                         errorCall, it,
                                                         shouldBe,
                                                         shouldSatisfy,
                                                         shouldThrow, shouldNotBe)
import Bio.NucleicAcid.Nucleotide (DNA(..), cNA)
import Bio.Tools.Sequence.CodonOptimization.Types (Organism(..))
import Bio.Protein.AminoAcid (AA(..))
import Debug.Trace (trace)
import Control.DeepSeq (force)

utilsSpec :: Spec
utilsSpec =
    describe "utilsSpec" $ do
        assembleSpec
        assembleRealSpec
        assembleWithGapSpec

        sliceSpec
        sliceOutOfBoundIndexSpec
        sliceWrongIndexesSpec

        buildOligSetSpec
        buildEmptyOligSetSpec
        buildOligSetForEmptySplittingSpec
        buildOligSetSplittingSpec
        buildOligSetWithGapSplittingSpec
        buildOligSetWithIncorrectSplittingSpec
        buildOligSetWithOutOfBoundSplittingSpec
        buildOligSetWithOneSplittingCoordinateSpec

assembleSpec :: Spec
assembleSpec =
    describe "assembleSpec" $
    it "assembleSpec" $ do
        let fwd = [Olig "AA" 0 2, Olig "TT" 2 4, Olig "GG" 4 6, Olig "CC" 6 8, Olig "AA" 8 10]
        let rsd = [Olig "AT" 1 3, Olig "CA" 3 5, Olig "GC" 5 7, Olig "TG" 7 9, Olig "AT" 9 11]
        let oligs = OligSet fwd rsd (OligSplitting [] [])
        assemble oligs `shouldBe` "AATTGGCCAAT"

assembleRealSpec :: Spec
assembleRealSpec =
    describe "assembleRealSpec'" $
    it "assembleRealSpec" $ do
        let coords = OligSplitting [(0, 60), (60, 120), (120, 180), (180, 240)] [(30, 90), (90, 150), (150, 210), (210, 270)]
        let oligs =
                OligSet
                    [ Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60
                    , Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120
                    , Olig "TGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGC" 120 180
                    , Olig "GGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACC" 180 240
                    ]
                    [ Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90
                    , Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150
                    , Olig "GGTCACCACGCTGCTCAGGCTGTACAGGCCGCTGCTCTGCAGCACGGCAGGGAAGGTGTG" 150 210
                    , Olig "GCTAGGCTTGTGGTTCACGTTGCAGATGTAGGTCTGGGTGCCCAGGCTGCTGCTAGGCAC" 210 270
                    ]
                    coords
        assemble oligs `shouldBe` "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGCGGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGCTGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGCGGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACCTACATCTGCAACGTGAACCACAAGCCTAGC"

assembleWithGapSpec :: Spec
assembleWithGapSpec =
    describe "assembleWithGapSpec" $
    it "assembleWithGapSpec" $ do
        let fwd = [Olig "AATT" 0 4, Olig "CCAA" 6 10]
        let rsd = [Olig "GCCA" 3 7, Olig "ATT" 8 11]
        let oligs = OligSet fwd rsd (OligSplitting [] [])
        assemble oligs `shouldBe` "AATTGGCCAAT"

sliceSpec :: Spec
sliceSpec =
    describe "sliceSpec" $
    it "should correct get slice from sequence" $ do
        let sequence = [0..19] :: [Integer]
        let res = slice 3 6 sequence
        res `shouldBe` [3, 4, 5]

sliceWrongIndexesSpec :: Spec
sliceWrongIndexesSpec =
    describe "sliceWrongIndexesSpec" $
    it "should return empty slice for wrong indexes" $ do
        let sequence = [0..19] :: [Integer]
        evaluate (slice 6 3 sequence)    `shouldThrow` errorCall "incorrect coordinates"
        evaluate (slice (-6) 0 sequence) `shouldThrow` errorCall "incorrect coordinates"
        slice 20 22 sequence `shouldBe` []
        slice 0 0 sequence `shouldBe` []
        evaluate (slice (-3) 2 sequence) `shouldThrow` errorCall "incorrect coordinates"

sliceOutOfBoundIndexSpec :: Spec
sliceOutOfBoundIndexSpec =
    describe "sliceOutOfBoundIndexSpec" $
    it "should return tail for out of bound indexes" $ do
        let sequence = [0..19] :: [Integer]
        slice 17 33 sequence `shouldBe` [17, 18, 19]
        slice 19 20 sequence `shouldBe` [19]

buildEmptyOligSetSpec :: Spec
buildEmptyOligSetSpec =
    describe "buildEmptyOligSetSpec" $
    it "should build empty oligs from splitting and sequence" $ do
        let splitting = OligSplitting [] []
        let dna = ""
        let res = buildOligSet splitting dna
        res `shouldBe` (OligSet [] [] splitting)

buildOligSetForEmptySplittingSpec :: Spec
buildOligSetForEmptySplittingSpec =
    describe "buildOligSetForEmptySplittingSpec" $
    it "should build empty oligs from empty splitting and sequence" $ do
        let splitting = OligSplitting [] []
        let dna = "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGGTCGACCGGCATGGCTTCCATGTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGTGGCAGGGGG"
        let res = buildOligSet splitting dna
        res `shouldBe` (OligSet [] [] splitting)

buildOligSetSplittingSpec :: Spec
buildOligSetSplittingSpec =
    describe "buildOligSetSplittingSpec" $
    it "should correct build oligs from splitting and sequence" $ do
        let splitting = OligSplitting [(0, 57), (57, 114)] [(29, 86), (86, 123)]
        let dna = "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGGTCGACCGGCATGGCTTCCATGTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGTGGCAGGGGG"
        let res = buildOligSet splitting dna
        res `shouldBe` OligSet [ Olig "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGGTCGACC" 0 57
                               , Olig "GGCATGGCTTCCATGTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGT" 57 114]
                               [ Olig (reverse $ translate "GCTGCTGCTGTGGGTGCCTGGGTCGACCGGCATGGCTTCCATGTCGGCAGAATGCTT") 29 86
                               , Olig (reverse $ translate "AATGAATTACAACAGTACTGCGATGAGTGGCAGGGGG") 86 123]
                                splitting
        assemble res `shouldBe` dna

buildOligSetSpec :: Spec
buildOligSetSpec =
    describe "buildOligSetSpec" $
    it "should correct build oligs from splitting and sequence" $ do
        let splitting = OligSplitting [(0, 56), (56, 112), (112, 168), (168, 224), (224, 280)]
                                                               [(28, 84), (84, 140), (140, 196), (196, 252), (252, 306)]
        let dna = "GACATACAAATGACTCAAAGTCCATCTACTCTATCCGCGAGTGTCGGCGACCGCGTAACTATTACGTGCAGGGCTTCACAAAGCATCGGTTCGGCTTTAGCATGGTATCAGCAGAAGCCTGGGAAAGCTCCTAAGTTACTGATCTATAAGGCAAGTGCCCTGGAGAACGGTGTTCCGTCTAGGTTTTCGGGCTCTGGTAGTGGGACCGAGTTCACACTGACAATAAGCAGTCTCCAACCCGATGATTTCGCCACCTACTACTGCCAGCACCTGACCTTCGGACAAGGGACGAGGTTGGAAATCAAA"
        let res = buildOligSet splitting dna
        res `shouldBe` OligSet
                          [Olig "GACATACAAATGACTCAAAGTCCATCTACTCTATCCGCGAGTGTCGGCGACCGCGT" 0 56,
                           Olig "AACTATTACGTGCAGGGCTTCACAAAGCATCGGTTCGGCTTTAGCATGGTATCAGC" 56 112,
                           Olig "AGAAGCCTGGGAAAGCTCCTAAGTTACTGATCTATAAGGCAAGTGCCCTGGAGAAC" 112 168,
                           Olig "GGTGTTCCGTCTAGGTTTTCGGGCTCTGGTAGTGGGACCGAGTTCACACTGACAAT" 168 224,
                           Olig "AAGCAGTCTCCAACCCGATGATTTCGCCACCTACTACTGCCAGCACCTGACCTTCG" 224 280]
                          [Olig "GCTTTGTGAAGCCCTGCACGTAATAGTTACGCGGTCGCCGACACTCGCGGATAGAG" 28 84,
                           Olig "AGTAACTTAGGAGCTTTCCCAGGCTTCTGCTGATACCATGCTAAAGCCGAACCGAT" 84 140,
                           Olig "CAGAGCCCGAAAACCTAGACGGAACACCGTTCTCCAGGGCACTTGCCTTATAGATC" 140 196,
                           Olig "GGCGAAATCATCGGGTTGGAGACTGCTTATTGTCAGTGTGAACTCGGTCCCACTAC" 196 252,
                           Olig "TTTGATTTCCAACCTCGTCCCTTGTCCGAAGGTCAGGTGCTGGCAGTAGTAGGT" 252 306]
                         splitting
        assemble res `shouldBe` dna

buildOligSetWithGapSplittingSpec :: Spec
buildOligSetWithGapSplittingSpec =
    describe "buildOligSetWithGapSplittingSpec" $
    it "should build oligs for splitting with gap and sequence" $ do
        let splitting = OligSplitting [(0, 57), (60, 117)] [(26, 83), (86, 123)]
        let dna = "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGGTCGACCGGCATGGCTTCCATGTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGTGGCAGGGGG"
        let res = buildOligSet splitting dna
        res `shouldBe` OligSet [ Olig "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGGTCGACC" 0 57
                               , Olig "ATGGCTTCCATGTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGTGGC" 60 117]
                               [ Olig (reverse $ translate "GGTGCTGCTGCTGTGGGTGCCTGGGTCGACCGGCATGGCTTCCATGTCGGCAGAATG") 26 83
                               , Olig (reverse $ translate "AATGAATTACAACAGTACTGCGATGAGTGGCAGGGGG") 86 123]
                                splitting
        assemble res `shouldBe` dna

buildOligSetWithOutOfBoundSplittingSpec :: Spec
buildOligSetWithOutOfBoundSplittingSpec =
    describe "buildOligSetWithOutOfBoundSplittingSpec" $
    it "should build oligs for splitting with out of bound end coordinate and sequence" $ do
        let splitting = OligSplitting [(0, 57), (60, 117)] [(26, 83), (86, 224)]
        let dna = "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGGTCGACCGGCATGGCTTCCATGTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGTGGCAGGGGG"
        let res = buildOligSet splitting dna
        res `shouldBe` OligSet [ Olig "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGGTCGACC" 0 57
                               , Olig "ATGGCTTCCATGTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGTGGC" 60 117]
                               [ Olig (reverse $ translate "GGTGCTGCTGCTGTGGGTGCCTGGGTCGACCGGCATGGCTTCCATGTCGGCAGAATG") 26 83
                               , Olig (reverse $ translate "AATGAATTACAACAGTACTGCGATGAGTGGCAGGGGG") 86 224]
                                splitting
        assemble res `shouldBe` dna

buildOligSetWithOneSplittingCoordinateSpec :: Spec
buildOligSetWithOneSplittingCoordinateSpec =
    describe "buildOligSetWithOneSplittingCoordinateSpec'" $
    it "should build oligs for splitting with one coordinate and sequence" $ do
        let splitting = OligSplitting [(0, 123)] [(26, 123)]
        let dna = "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGGTCGACCGGCATGGCTTCCATGTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGTGGCAGGGGG"
        let res = buildOligSet splitting dna
        res `shouldBe` OligSet [ Olig dna 0 123]
                               [ Olig (reverse $ translate $ drop 26 dna) 26 123]
                                splitting
        assemble res `shouldBe` dna

buildOligSetWithIncorrectSplittingSpec :: Spec
buildOligSetWithIncorrectSplittingSpec =
    describe "buildOligSetWithIncorrectSplittingSpec" $
    it "should build oligs for incorrect splitting and sequence" $ do
        let splitting = OligSplitting [(-10, 569)] [(29, 670)]
        let dna = "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGGTCGACCGGCATGGCTTCCATGTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGTGGCAGGGGG"
        let res = evaluate (buildOligSet splitting dna) `shouldThrow` errorCall "incorrect coordinates"

        let splitting' = OligSplitting [(10, 5)] [(29, 670)]
        (evaluate . force) (buildOligSet splitting dna) `shouldThrow` errorCall "incorrect coordinates"