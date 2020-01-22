module OligoDesigner.SpecOligoDesignerOptimizer where

import Test.Hspec (Spec, shouldBe, it, describe, shouldThrow, errorCall, shouldSatisfy)
import Bio.Tools.Sequence.OligoDesigner.Scorer (commonScore)
import Bio.Tools.Sequence.OligoDesigner.Optimizer (maxPairMutationIndexes, minPairMutationIndexes, mutationIndexes,
                                                   mutateSlice, mutate, minMaxOptimize)
import Bio.Tools.Sequence.OligoDesigner.Types     (Olig(..), MatrixCell(..), OligBounds, OligSplitting(..), OligSet(..),
                                                    OligSplittingConfig(..), OligoDesignerConfig(..))
import Data.Matrix (matrix)
import System.Random (mkStdGen)
import Control.Monad.State (evalState)
import Bio.Tools.Sequence.CodonOptimization.Types (Organism(..))
import Data.Default (def)
import Control.Exception (evaluate)
import Bio.Tools.Sequence.OligoDesigner.Utils (assemble, translate)

optimizerSpec :: Spec
optimizerSpec =
    describe "optimizerSpec" $ do
        minPairMutationIndexesSpec
        minPairMutationIndexesWithEvenIndexesSpec
        minPairMutationIndexesWithoutIntersectionSpec
        minPairMutationIndexesWithOneNKIntersectionSpec
        minPairMutationIndexesReverseIntersectionSpec
        minPairMutationIndexesWithBorderIndexesSpec

        maxPairMutationIndexesSpec
        maxPairMutationIndexesWithBorderIndexesSpec
        maxPairMutationIndexesWithEvenIndexesSpec

        mutationIndexesSpec
        mutationIndexesWithoutMinSpec

        mutateSliceSpec
        mutateSliceWhenThereIsNoVariantsSpec
        mutateSliceRealRandomSpec

        mutateSpec
        mutateFromStartSpec
        mutateInvalidIntervalSpec
        mutateOneAASpec

        minMaxOptimizeSpec

minPairMutationIndexesSpec :: Spec
minPairMutationIndexesSpec =
    describe "minPairMutationIndexes" $
    it "" $ do
        let matrixCell = MatrixCell (Olig "GCGGATTCAGCAGTT" 0 15) (Olig "ATTGGCGCATGCTTT" 10 25) (-34.993)
        let res = minPairMutationIndexes matrixCell
        res `shouldBe` [(4, 5)]

minPairMutationIndexesWithEvenIndexesSpec :: Spec
minPairMutationIndexesWithEvenIndexesSpec =
    describe "minPairMutationIndexesWithEvenIndexesSpec" $
    it "" $ do
        let matrixCell = MatrixCell (Olig "GCGGATTC" 0 8) (Olig "ATTGGCGC" 6 15) (-34.993)
        let res = minPairMutationIndexes matrixCell
        res `shouldBe` [(3, 3)]

minPairMutationIndexesWithBorderIndexesSpec :: Spec
minPairMutationIndexesWithBorderIndexesSpec =
    describe "minPairMutationIndexesWithEvenIndexesSpec" $
    it "" $ do
        let matrixCell = MatrixCell (Olig "GCGGATTCA" 0 9) (Olig "ATTGGCGCA" 6 15) (-34.993)
        let res = minPairMutationIndexes matrixCell
        res `shouldBe` [(3, 3)]

minPairMutationIndexesWithoutIntersectionSpec :: Spec
minPairMutationIndexesWithoutIntersectionSpec =
    describe "minPairMutationIndexesWithoutIntersectionSpec" $
    it "" $ do
        let matrixCell = MatrixCell (Olig "GCGGATTC" 0 8) (Olig "ATTGGCGC" 9 18) (-34.993)
        let res = minPairMutationIndexes matrixCell
        res `shouldBe` []

minPairMutationIndexesWithOneNKIntersectionSpec :: Spec
minPairMutationIndexesWithOneNKIntersectionSpec =
    describe "minPairMutationIndexesWithOneNKIntersectionSpec" $
    it "" $ do
        let matrixCell = MatrixCell (Olig "GCGGATTC" 0 8) (Olig "ATTGGCGC" 8 17) (-34.993)
        let res = minPairMutationIndexes matrixCell
        res `shouldBe` []

minPairMutationIndexesReverseIntersectionSpec :: Spec
minPairMutationIndexesReverseIntersectionSpec =
    describe "minPairMutationIndexesReverseIntersectionSpec" $
    it "" $ do
        let matrixCell = MatrixCell (Olig "GCGGATTC" 12 20) (Olig "ATTGGCGC" 1 10) (-34.993)
        let res = minPairMutationIndexes matrixCell
        res `shouldBe` []

maxPairMutationIndexesSpec :: Spec
maxPairMutationIndexesSpec =
    describe "maxPairMutationIndexesSpec" $
    it "" $ do
        let matrixCell = MatrixCell (Olig "GCGGATTC" 0 8) (Olig "ATTGGCGC" 8 17) (-34.993)
        let res = maxPairMutationIndexes matrixCell
        res `shouldBe` [(1, 3), (3, 6)]

maxPairMutationIndexesWithBorderIndexesSpec :: Spec
maxPairMutationIndexesWithBorderIndexesSpec =
    describe "maxPairMutationIndexesWithBorderIndexesSpec" $
    it "" $ do
        let matrixCell = MatrixCell (Olig "GCGGATTCA" 0 9) (Olig "ATTGGCAAA" 9 15) (-34.993)
        let res = maxPairMutationIndexes matrixCell
        res `shouldBe` [(1, 3), (4, 5)]

maxPairMutationIndexesWithEvenIndexesSpec :: Spec
maxPairMutationIndexesWithEvenIndexesSpec =
    describe "maxPairMutationIndexesWithEvenIndexesSpec" $
    it "" $ do
        let matrixCell = MatrixCell (Olig "GCGGATTC" 0 8) (Olig "ATTGGCAAATTT" 3 15) (-34.993)
        let res = maxPairMutationIndexes matrixCell
        res `shouldBe` [(1, 3), (2, 5)]

mutationIndexesSpec :: Spec
mutationIndexesSpec =
    describe "mutationIndexesSpec" $
    it "" $ do
        let mtx = matrix 5 5 generator
        let res = mutationIndexes mtx
        res `shouldBe` [(7, 9), (10, 17), (18, 30)]
  where
    generator :: (Int, Int) -> MatrixCell
    generator (1, 2) = MatrixCell (Olig "" 0 25) (Olig "" 20 36) (fromIntegral (-100500))
    generator (2, 5) = MatrixCell (Olig "" 29 51) (Olig "" 53 88) 100500
    generator (x, y) | abs (x - y) == 1 = MatrixCell (Olig "" x y) (Olig "" x y) (fromIntegral $ (-1) * x)
                     | otherwise        = MatrixCell (Olig "" x y) (Olig "" x y) (fromIntegral x)

mutationIndexesWithoutMinSpec :: Spec
mutationIndexesWithoutMinSpec =
    describe "mutationIndexesWithoutMinSpec" $
    it "" $ do
        let mtx = matrix 5 5 generator
        let res = mutationIndexes mtx
        res `shouldBe` [(10, 17), (18, 30)]
  where
    generator :: (Int, Int) -> MatrixCell
    generator (1, 2) = MatrixCell (Olig "" 0 25) (Olig "" 35 51) (fromIntegral (-100500))
    generator (2, 5) = MatrixCell (Olig "" 29 51) (Olig "" 53 88) 100500
    generator (x, y) | abs (x - y) == 1 = MatrixCell (Olig "" x y) (Olig "" x y) (fromIntegral $ (-1) * x)
                     | otherwise        = MatrixCell (Olig "" x y) (Olig "" x y) (fromIntegral x)

mutateSliceSpec :: Spec
mutateSliceSpec =
    describe "mutateSliceSpec" $
    it "" $ do
        let gen = mkStdGen 4
        let res = evalState (mutateSlice CHO "AATATGCAT") gen
        res `shouldBe` ["AATATGCAT", "AACATGCAT", "AATATGCAC"]

mutateSliceWhenThereIsNoVariantsSpec :: Spec
mutateSliceWhenThereIsNoVariantsSpec =
    describe "mutateSliceWhenThereIsNoVariantsSpec" $
    it "" $ do
        let gen = mkStdGen 4
        let res = evalState (mutateSlice CHO "ATGTGG") gen
        res `shouldBe` ["ATGTGG"]

mutateSliceRealRandomSpec :: Spec
mutateSliceRealRandomSpec =
    describe "mutateSliceRealRandomSpec" $
    it "" $ do
        let gen = mkStdGen 9
        let res = evalState (mutateSlice CHO "TCTTTGCCGAACGAGGGCATG") gen
        res `shouldBe` ["TCTTTGCCGAACGAGGGCATG", "AGCTTGCCGAACGAGGGCATG", "TCTCTCCCGAACGAGGGCATG", "TCTTTGCCAAACGAGGGCATG",
            "TCTTTGCCGAATGAGGGCATG", "TCTTTGCCGAACGAAGGCATG", "TCTTTGCCGAACGAGGGAATG"]

mutateSpec :: Spec
mutateSpec =
    describe "mutateSpec" $
    it "" $ do
        let dna = "ATGGAGACC" ++ "AATATGCAT" ++ "GACACCCTGCTGCTGTGGGTGCTGCTGCTG"
        let gen = mkStdGen 4
        let res = evalState (mutate def dna (4, 6)) gen
        res `shouldBe` ["ATGGAGACCAATATGCATGACACCCTGCTGCTGTGGGTGCTGCTGCTG",
                        "ATGGAGACCAACATGCATGACACCCTGCTGCTGTGGGTGCTGCTGCTG",
                        "ATGGAGACCAATATGCACGACACCCTGCTGCTGTGGGTGCTGCTGCTG"]

mutateFromStartSpec :: Spec
mutateFromStartSpec =
    describe "mutateFromStartSpec" $
    it "" $ do
        let dna = "AATATGCAT" ++ "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTG"
        let gen = mkStdGen 4
        let res = evalState (mutate def dna (1, 3)) gen
        res `shouldBe` ["AATATGCATATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTG",
                        "AACATGCATATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTG",
                        "AATATGCACATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTG"]

mutateOneAASpec :: Spec
mutateOneAASpec =
    describe "mutateOneAASpec" $
    it "" $ do
        let dna = "ATGGAGACC" ++ "AAT" ++ "ATGCATGACACCCTGCTGCTGTGGGTGCTGCTGCTG"
        let gen = mkStdGen 4
        let res = evalState (mutate def dna (4, 4)) gen
        res `shouldBe` ["ATGGAGACCAATATGCATGACACCCTGCTGCTGTGGGTGCTGCTGCTG",
                        "ATGGAGACCAACATGCATGACACCCTGCTGCTGTGGGTGCTGCTGCTG"]

mutateInvalidIntervalSpec :: Spec
mutateInvalidIntervalSpec =
    describe "mutateInvalidIntervalSpec" $
    it "" $ do
        let gen = mkStdGen 4
        evaluate (evalState (mutate def "AATATGCATATG" (3, 1)) gen) `shouldThrow` errorCall "invalid interval for mutation: (3,1)"
        evaluate (evalState (mutate def "AATATGCATATG" (-1, 1)) gen) `shouldThrow` errorCall "invalid interval for mutation: (-1,1)"
        evaluate (evalState (mutate def "AATATGCATATG" (3, -10)) gen) `shouldThrow` errorCall "invalid interval for mutation: (3,-10)"
        evaluate (evalState (mutate def "AATATGCATATG" (3, 5)) gen) `shouldThrow` errorCall "invalid interval for mutation: (3,5)"
        evaluate (evalState (mutate def "AATATGCATATG" (5, 40)) gen) `shouldThrow` errorCall "invalid interval for mutation: (5,40)"

minMaxOptimizeSpec :: Spec
minMaxOptimizeSpec =
    describe "minMaxOptimizeScore" $
    it "" $ do
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
        let conf = OligoDesignerConfig def 0 (OligSplittingConfig 60 1 10)
        let gen = mkStdGen 499
        let res = evalState (minMaxOptimize conf oligs) gen

        assemble res `shouldBe` "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCTCTAAGAGCACCAGCGGCGGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGCTGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGCGGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACCTACATCTGCAACGTGAACCACAAGCCTAGC"
        res `shouldBe` OligSet
                           [ Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCTCTAAGAGCACCAGCGGC" 0 60
                           , Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120
                           , Olig "TGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGC" 120 180
                           , Olig "GGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACC" 180 240
                           ]
                           [ Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTAGAGCTAGGGGCCAG" 30 90
                           , Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150
                           , Olig "GGTCACCACGCTGCTCAGGCTGTACAGGCCGCTGCTCTGCAGCACGGCAGGGAAGGTGTG" 150 210
                           , Olig "GCTAGGCTTGTGGTTCACGTTGCAGATGTAGGTCTGGGTGCCCAGGCTGCTGCTAGGCAC" 210 270
                           ]
                           coords
        commonScore conf res `shouldSatisfy` (>= commonScore conf oligs)