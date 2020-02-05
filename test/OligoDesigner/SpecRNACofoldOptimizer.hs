module OligoDesigner.SpecRNACofoldOptimizer where

import Test.Hspec (Spec, shouldBe, it, describe, shouldThrow, errorCall, shouldSatisfy)
import Bio.Tools.Sequence.OligoDesigner.Scorer (commonScore)
import Bio.Tools.Sequence.OligoDesigner.Optimizer.RNACofoldOptimizer (maxPairMutationIndexes, minPairMutationIndexes,
    mutationIndexes, rnaOptimize)
import Bio.Tools.Sequence.OligoDesigner.Types     (MatrixCell(..), OligBounds, OligSplitting(..), OligSet(..),
                                                    OligsSplittingConfig(..), OligsDesignerConfig(..), Olig(..), OligLight(..))
import Data.Matrix (matrix)
import System.Random (mkStdGen)
import Control.Monad.State (evalState)
import Bio.Tools.Sequence.CodonOptimization.Types (Organism(..))
import Data.Default (def)
import Control.Exception (evaluate)
import Bio.Tools.Sequence.OligoDesigner.Utils (assemble, translate)

rnaOptimizerSpec :: Spec
rnaOptimizerSpec =
    describe "rnaOptimizerSpec" $ do
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

        rnaOptimizeSpec
        
minPairMutationIndexesSpec :: Spec
minPairMutationIndexesSpec =
    describe "minPairMutationIndexes" $
    it "" $ do
        let matrixCell = MatrixCell (OligLight "GCGGATTCAGCAGTT" (Olig "GCGGATTCAGCAGTT" 0 15)) 
                                    (OligLight "ATTGGCGCATGCTTT" (Olig "ATTGGCGCATGCTTT" 10 25)) (-34.993)
        let res = minPairMutationIndexes matrixCell
        res `shouldBe` [(4, 5)]

minPairMutationIndexesWithEvenIndexesSpec :: Spec
minPairMutationIndexesWithEvenIndexesSpec =
    describe "minPairMutationIndexesWithEvenIndexesSpec" $
    it "" $ do
        let matrixCell = MatrixCell (OligLight "" (Olig "" 0 8)) (OligLight "" (Olig "" 6 15)) (-34.993)
        let res = minPairMutationIndexes matrixCell
        res `shouldBe` [(3, 3)]

minPairMutationIndexesWithBorderIndexesSpec :: Spec
minPairMutationIndexesWithBorderIndexesSpec =
    describe "minPairMutationIndexesWithEvenIndexesSpec" $
    it "" $ do
        let matrixCell = MatrixCell (OligLight "" (Olig "" 0 9)) (OligLight "" (Olig "" 6 15)) (-34.993)
        let res = minPairMutationIndexes matrixCell
        res `shouldBe` [(3, 3)]

minPairMutationIndexesWithoutIntersectionSpec :: Spec
minPairMutationIndexesWithoutIntersectionSpec =
    describe "minPairMutationIndexesWithoutIntersectionSpec" $
    it "" $ do
        let matrixCell = MatrixCell (OligLight "" (Olig "" 0 8)) (OligLight "" (Olig "" 9 18)) (-34.993)
        let res = minPairMutationIndexes matrixCell
        res `shouldBe` []

minPairMutationIndexesWithOneNKIntersectionSpec :: Spec
minPairMutationIndexesWithOneNKIntersectionSpec =
    describe "minPairMutationIndexesWithOneNKIntersectionSpec" $
    it "" $ do
        let matrixCell = MatrixCell (OligLight "" (Olig "" 0 8)) (OligLight "" (Olig "" 8 17)) (-34.993)
        let res = minPairMutationIndexes matrixCell
        res `shouldBe` []

minPairMutationIndexesReverseIntersectionSpec :: Spec
minPairMutationIndexesReverseIntersectionSpec =
    describe "minPairMutationIndexesReverseIntersectionSpec" $
    it "" $ do
        let matrixCell = MatrixCell (OligLight "" (Olig "" 12 20)) (OligLight "" (Olig "" 1 10)) (-34.993)
        let res = minPairMutationIndexes matrixCell
        res `shouldBe` []

maxPairMutationIndexesSpec :: Spec
maxPairMutationIndexesSpec =
    describe "maxPairMutationIndexesSpec" $
    it "" $ do
        let matrixCell = MatrixCell (OligLight "" (Olig "" 0 8)) (OligLight "" (Olig "" 8 17)) (-34.993)
        let res = maxPairMutationIndexes matrixCell
        res `shouldBe` [(1, 3), (3, 6)]

maxPairMutationIndexesWithBorderIndexesSpec :: Spec
maxPairMutationIndexesWithBorderIndexesSpec =
    describe "maxPairMutationIndexesWithBorderIndexesSpec" $
    it "" $ do
        let matrixCell = MatrixCell (OligLight "" (Olig "" 0 9)) (OligLight "" (Olig "" 9 15)) (-34.993)
        let res = maxPairMutationIndexes matrixCell
        res `shouldBe` [(1, 3), (4, 5)]

maxPairMutationIndexesWithEvenIndexesSpec :: Spec
maxPairMutationIndexesWithEvenIndexesSpec =
    describe "maxPairMutationIndexesWithEvenIndexesSpec" $
    it "" $ do
        let matrixCell = MatrixCell (OligLight "" (Olig "" 0 8)) (OligLight "" (Olig "" 3 15)) (-34.993)
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
    generator (1, 2) = MatrixCell (OligLight "" (Olig "" 0 25)) (OligLight "" (Olig "" 20 36)) (fromIntegral (-100500))
    generator (2, 5) = MatrixCell (OligLight "" (Olig "" 29 51)) (OligLight "" (Olig "" 53 88)) 100500
    generator (x, y) | abs (x - y) == 1 = MatrixCell (OligLight "" (Olig "" x y)) (OligLight "" (Olig "" x y)) (fromIntegral $ (-1) * x)
                     | otherwise        = MatrixCell (OligLight "" (Olig "" x y)) (OligLight "" (Olig "" x y)) (fromIntegral x)

mutationIndexesWithoutMinSpec :: Spec
mutationIndexesWithoutMinSpec =
    describe "mutationIndexesWithoutMinSpec" $
    it "" $ do
        let mtx = matrix 5 5 generator
        let res = mutationIndexes mtx
        res `shouldBe` [(10, 17), (18, 30)]
  where
    generator :: (Int, Int) -> MatrixCell
    generator (1, 2) = MatrixCell (OligLight "" (Olig "" 0 25)) (OligLight "" (Olig "" 35 51)) (fromIntegral (-100500))
    generator (2, 5) = MatrixCell (OligLight "" (Olig "" 29 51)) (OligLight "" (Olig "" 53 88)) 100500
    generator (x, y) | abs (x - y) == 1 = MatrixCell (OligLight "" (Olig "" x y)) (OligLight "" (Olig "" x y)) (fromIntegral $ (-1) * x)
                     | otherwise        = MatrixCell (OligLight "" (Olig "" x y)) (OligLight "" (Olig "" x y)) (fromIntegral x)

rnaOptimizeSpec :: Spec
rnaOptimizeSpec =
    describe "rnaOptimizeSpec" $
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
        let conf = OligsDesignerConfig def (OligsSplittingConfig 60 1 10) 0 0 0 0
        let gen = mkStdGen 499
        let res = evalState (rnaOptimize conf oligs) gen

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