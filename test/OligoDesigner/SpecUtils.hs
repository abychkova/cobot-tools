{-# LANGUAGE OverloadedStrings #-}

module OligoDesigner.SpecUtils where

import           Bio.Tools.Sequence.OligoDesigner.Types (Olig (..),
                                                         OligSet (..),
                                                         OligSplitting (..))
import           Bio.Tools.Sequence.OligoDesigner.Utils (assemble,
                                                         weightedRandom,
                                                         slice,
                                                         buildOligSet,
                                                         randomCodon, oneMutation, translate)
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

        weightedRandomSimpleSpec
        weightedRandomMultipleCallSpec
        weightedRandomSpec
        weightedRandomAllZeroSpec
        weightedRandomForUnsortedSpec
        weightedRandomForEmptySpec

        sliceSpec
        sliceOutOfBoundIndexSpec
        sliceWrongIndexesSpec

        buildEmptyOligSetSpec
        buildOligSetForEmptySplittingSpec
        buildOligSetSplittingSpec
        buildOligSetWithGapSplittingSpec
        buildOligSetWithIncorrectSplittingSpec
        buildOligSetWithOutOfBoundSplittingSpec
        buildOligSetWithOneSplittingCoordinateSpec

        randomCodonSpec
        randomCodonForAAWithOneCodonSpec
        randomCodonWeightedSpec
        randomCodonForDifferentOrganismSpec

        oneMutationSpec
        oneMutationForAAWithOneCodonSpec

assembleSpec :: Spec
assembleSpec =
    describe "assembleSpec" $
    it "assembleSpec" $ do
        let fwd = [Olig "AA" 0 2, Olig "TT" 2 4, Olig "GG" 4 6, Olig "CC" 6 8, Olig "AA" 8 10]
        let rsd = [Olig "TA" 1 3, Olig "AC" 3 5, Olig "CG" 5 7, Olig "GT" 7 9, Olig "TA" 9 11]
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
                    [ Olig "GACCGGGGATCGTCGTTCTCGTGGTCGCCGCCGTGGCGGCGGGACCCGACGGACCACTTC" 30 90
                    , Olig "CTGATGAAGGGACTCGGACACTGGCACTCGACCTTGTCGCCGCGGGACTGGTCGCCGCAC" 90 150
                    , Olig "GTGTGGAAGGGACGGCACGACGTCTCGTCGCCGGACATGTCGGACTCGTCGCACCACTGG" 150 210
                    , Olig "CACGGATCGTCGTCGGACCCGTGGGTCTGGATGTAGACGTTGCACTTGGTGTTCGGATCG" 210 270
                    ]
                    coords
        assemble oligs `shouldBe` "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGCGGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGCTGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGCGGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACCTACATCTGCAACGTGAACCACAAGCCTAGC"

assembleWithGapSpec :: Spec
assembleWithGapSpec =
    describe "assembleWithGapSpec" $
    it "assembleWithGapSpec" $ do
        let fwd = [Olig "AATT" 0 4, Olig "CCAA" 6 10]
        let rsd = [Olig "ACCG" 3 7, Olig "TTA" 8 11]
        let oligs = OligSet fwd rsd (OligSplitting [] [])
        assemble oligs `shouldBe` "AATTGGCCAAT"

weightedRandomSimpleSpec :: Spec
weightedRandomSimpleSpec =
    describe "weightedRandomSimpleSpec" $
    it "should return weigthed random element" $ do
        gen <- getStdGen
        let arg = [(1, 1), (2, 0)] :: [(Integer, Double)]
        let res = evalState (weightedRandom arg) gen
        res `shouldBe` 1
        let arg2 = [(1, 0.99999999), (2, 0.00000001)] :: [(Integer, Double)]
        let res2 = evalState (weightedRandom arg2) gen
        res2 `shouldBe` 1
        let arg3 = [(1, 0), (2, 1)] :: [(Integer, Double)]
        let res3 = evalState (weightedRandom arg3) gen
        res3 `shouldBe` 2

weightedRandomMultipleCallSpec :: Spec
weightedRandomMultipleCallSpec =
    describe "weightedRandomMultipleCallSpec" $
    it "random elements ahould be different for each call" $ do
        gen <- getStdGen
        let arg = [(1, 1), (2, 1)] :: [(Integer, Double)]
        let res = runWeightedRandomNTimes gen arg 3
        length res `shouldSatisfy` (> (length $ nub res))

weightedRandomSpec :: Spec
weightedRandomSpec =
    describe "weightedRandomSpec" $
    it "should return weigthed random element" $ do
        gen <- getStdGen
        let res = runWeightedRandomNTimes gen [(3, 7), (4, 300), (1, 1600)] 100
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
        gen <- getStdGen
        let res = runWeightedRandomNTimes gen [(4, 300), (3, 7), (1, 1600)] 100
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
        gen <- getStdGen
        let arg = [(1, 0), (2, 0)] :: [(Integer, Double)]
        let res = evalState (weightedRandom arg) gen
        res `shouldBe` 1

weightedRandomForEmptySpec :: Spec
weightedRandomForEmptySpec =
    describe "weightedRandomForEmptySpec" $
    it "should return error for empty argument" $ do
        let empty = [] :: [(Double, Double)]
        gen <- getStdGen
        evaluate (evalState (weightedRandom empty) gen) `shouldThrow` errorCall "cannot get random for empty array"

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
                               [ Olig (translate "GCTGCTGCTGTGGGTGCCTGGGTCGACCGGCATGGCTTCCATGTCGGCAGAATGCTT") 29 86
                               , Olig (translate "AATGAATTACAACAGTACTGCGATGAGTGGCAGGGGG") 86 123]
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
                               [ Olig (translate "GGTGCTGCTGCTGTGGGTGCCTGGGTCGACCGGCATGGCTTCCATGTCGGCAGAATG") 26 83
                               , Olig (translate "AATGAATTACAACAGTACTGCGATGAGTGGCAGGGGG") 86 123]
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
                               [ Olig (translate "GGTGCTGCTGCTGTGGGTGCCTGGGTCGACCGGCATGGCTTCCATGTCGGCAGAATG") 26 83
                               , Olig (translate "AATGAATTACAACAGTACTGCGATGAGTGGCAGGGGG") 86 224]
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
                               [ Olig (translate $ drop 26 dna) 26 123]
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

randomCodonSpec :: Spec
randomCodonSpec =
    describe "randomCodonSpec" $
    it "should return two different codons for HIS for defferent random" $ do
        let gen = mkStdGen 3
        let res = evalState (randomCodon CHO HIS) gen
        let res2 = evalState (randomCodon CHO HIS) gen
        res `shouldBe` [DC, DA, DT]
        res `shouldBe` res2

        let gen' = mkStdGen 5
        let res' = evalState (randomCodon CHO HIS) gen'
        res' `shouldBe` [DC, DA, DC]

randomCodonForAAWithOneCodonSpec :: Spec
randomCodonForAAWithOneCodonSpec =
    describe "randomCodonForAAWithOneCodonSpec" $
    it "should return only codon for MET for defferent random" $ do
        let gen = mkStdGen 3
        let res = evalState (randomCodon CHO MET) gen
        let gen' = mkStdGen 500
        let res' = evalState (randomCodon CHO MET) gen'
        res `shouldBe` [DA, DT, DG]
        res' `shouldBe` [DA, DT, DG]

randomCodonWeightedSpec :: Spec
randomCodonWeightedSpec =
    describe "randomCodonWeightedSpec" $
    it "should return most weghted codon for LEU maximum timmes" $ do
        let gens = [mkStdGen x | x <- [0,10..300]]
        let res = [evalState (randomCodon CHO LEU) gen | gen <- gens]
        let ctg = length (filter (== [DC, DT, DG]) res)
        let tta = length (filter (== [DT, DT, DA]) res)
        let ttg = length (filter (== [DT, DT, DG]) res)
        let ctt = length (filter (== [DC, DT, DT]) res)
        let ctc = length (filter (== [DC, DT, DC]) res)
        let cta = length (filter (== [DC, DT, DA]) res)

        ctg `shouldSatisfy` (>= tta)
        ctg `shouldSatisfy` (>= ttg)
        ctg `shouldSatisfy` (>= ctt)
        ctg `shouldSatisfy` (>= ctc)
        ctg `shouldSatisfy` (>= cta)

randomCodonForDifferentOrganismSpec :: Spec
randomCodonForDifferentOrganismSpec =
    describe "randomCodonForDifferentOrganismSpec" $
    it "should return defferent codon for ILE for different organism" $ do
        let gen = mkStdGen 7
        let humanRes = evalState (randomCodon Human ILE) gen
        let ecoliRes = evalState (randomCodon EColi ILE) gen
        humanRes `shouldNotBe` ecoliRes

oneMutationSpec :: Spec
oneMutationSpec =
    describe "oneMutationSpec" $
    it "should return another codon for LYS (aa with 2 codons) for any random value" $ do
        gen <- getStdGen
        let res = evalState (oneMutation CHO [DA, DA, DA]) gen
        res `shouldBe` [DA, DA, DG]

oneMutationForAAWithOneCodonSpec :: Spec
oneMutationForAAWithOneCodonSpec =
    describe "oneMutationForAAWithOneCodonSpec" $
    it "should return the same codon for TRP (aa with 1 codon) for any random value" $ do
        gen <- getStdGen
        let res = evalState (oneMutation CHO [DT, DG, DG]) gen
        res `shouldBe` [DT, DG, DG]

runWeightedRandomNTimes :: StdGen -> [(Integer, Double)] -> Int -> [Integer]
runWeightedRandomNTimes gen arg n = evalState (helper n []) gen
  where
    helper :: Int -> [Integer] -> State StdGen [Integer]
    helper 0 acc = return acc
    helper count acc = do
        random <- get
        let (res, newGen) = runState (weightedRandom arg) random
        put newGen
        helper (count - 1) (res : acc)