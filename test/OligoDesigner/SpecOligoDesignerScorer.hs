{-# LANGUAGE OverloadedStrings #-}

module OligoDesigner.SpecOligoDesignerScorer where

import           Bio.Tools.Sequence.OligoDesigner.Scorer (rnaScore, rnaMatrix, rnaMatrixScore)
import           Bio.Tools.Sequence.OligoDesigner.Types  (Olig (..),
                                                          OligSet (..),
                                                          OligSplitting (..), MatrixCell(..))
import           Test.Hspec                              (Spec, describe, it,
                                                          shouldBe)
import Data.Matrix (matrix)
import Bio.Tools.Sequence.OligoDesigner.Utils (assemble)

oligoDesignerScoreSpec :: Spec
oligoDesignerScoreSpec = describe "Oligo-Designer score spec" $ do
    scoreRNAFoldingSpec
    scoreRNAFoldingSpec2
    rnaMatrixSpec
    rnaMatrixSimpleSpec
    rnaMatrixScoreSpec


scoreRNAFoldingSpec :: Spec
scoreRNAFoldingSpec =
    describe "scoreRNAFoldingSpec" $
    it "should correct score oligs" $ do
        let coords = OligSplitting [(0, 60), (60, 120)] [(30, 90), (90, 150)]
        let oligs =
                OligSet
                    [ Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60
                    , Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120
                    ]
                    [ Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90
                    , Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150
                    ]
                    coords
        rnaScore oligs `shouldBe` -55.299995

scoreRNAFoldingSpec2 :: Spec
scoreRNAFoldingSpec2 =
    describe "scoreRNAFoldingSpec2" $
    it "should correct score oligs" $ do
         let coords = OligSplitting [(0, 60), (60, 120), (120, 180), (180, 240), (240, 300), (300, 360)]
                                    [(30, 90), (90, 150), (150, 210), (210, 270), (270, 330), (330, 390)]
         let oligs =
                 OligSet
                     [ Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60
                     , Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120
                     , Olig "TGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGC" 120 180
                     , Olig "GGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACC" 180 240
                     , Olig "TACATCTGCAACGTGAACCACAAGCCTAGCAACACCAAGGTGGACAAGAAGGTGGAGCCC" 240 300
                     , Olig "AAGAGCTGCGACCGTACGCACACCTGCCCCCCCTGCCCCGCCCCCGAGCTGCTGGGCGGC" 300 360
                     ]
                     [ Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90
                     , Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150
                     , Olig "GGTCACCACGCTGCTCAGGCTGTACAGGCCGCTGCTCTGCAGCACGGCAGGGAAGGTGTG" 150 210
                     , Olig "GCTAGGCTTGTGGTTCACGTTGCAGATGTAGGTCTGGGTGCCCAGGCTGCTGCTAGGCAC" 210 270
                     , Olig "GGGGCAGGTGTGCGTACGGTCGCAGCTCTTGGGCTCCACCTTCTTGTCCACCTTGGTGTT" 270 330
                     , Olig "AGGCTTAGGAGGGAACAGGAACACGCTAGGGCCGCCCAGCAGCTCGGGGGCGGGGCAGGG" 330 390
                     ]
                     coords
         rnaScore oligs `shouldBe` -71.1

rnaMatrixSpec :: Spec
rnaMatrixSpec =
    describe "rnaMatrixSpec" $
    it "" $ do
        let coords = OligSplitting [(0, 60), (60, 120)] [(30, 90), (90, 150)]
        let oligs =
                OligSet
                    [ Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60
                    , Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120
                    ]
                    [ Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90
                    , Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150
                    ]
                    coords
        let res = rnaMatrix oligs
        res `shouldBe` matrix 4 4 generatorRealMatrix

rnaMatrixSimpleSpec :: Spec
rnaMatrixSimpleSpec =
    describe "rnaMatrixSimpleSpec" $
    it "" $ do
        let coords = OligSplitting [(0, 6), (6, 12)] [(3, 9), (9, 15)]
        let oligs =
                OligSet
                    [ Olig "AATGGC" 0 6
                    , Olig "GACTGA" 6 12
                    ]
                    [ Olig "CCGCTG" 3 9
                    , Olig "ACTTTC" 9 15
                    ]
                    coords
        let res = rnaMatrix oligs
        res `shouldBe` matrix 4 4 generator
      where
        generator :: (Int, Int) -> MatrixCell
        generator (1, 1) = MatrixCell (Olig "AATGGC" 0 6)
                                      (Olig "AATGGC" 0 6) 0
        generator (1, 2) = MatrixCell (Olig "AATGGC" 0 6)
                                      (Olig "CCGCTG" 3 9) (-2)
        generator (1, 3) = MatrixCell (Olig "AATGGC" 0 6)
                                      (Olig "GACTGA" 6 12) 0
        generator (1, 4) = MatrixCell (Olig "AATGGC" 0 6)
                                      (Olig "ACTTTC" 9 15) 0
        generator (2, 2) = MatrixCell (Olig "CCGCTG" 3 9)
                                      (Olig "CCGCTG" 3 9) (-2.1)
        generator (2, 3) = MatrixCell (Olig "CCGCTG" 3 9)
                                      (Olig "GACTGA" 6 12) 0
        generator (2, 4) = MatrixCell (Olig "CCGCTG" 3 9)
                                      (Olig "ACTTTC" 9 15) 0
        generator (3, 3) = MatrixCell (Olig "GACTGA" 6 12)
                                      (Olig "GACTGA" 6 12) 0
        generator (3, 4) = MatrixCell (Olig "GACTGA" 6 12)
                                      (Olig "ACTTTC" 9 15) 0
        generator (4, 4) = MatrixCell (Olig "ACTTTC" 9 15)
                                      (Olig "ACTTTC" 9 15) 0
        generator (x, y) = MatrixCell (Olig "" 0 0 ) (Olig "" 0 0 ) 0

rnaMatrixScoreSpec :: Spec
rnaMatrixScoreSpec =
    describe "rnaMatrixScoreSpec" $
    it "" $ do
        let mtx = matrix 4 4 generatorRealMatrix
        let res = rnaMatrixScore mtx
        res `shouldBe` -55.299995

generatorRealMatrix :: (Int, Int) -> MatrixCell
generatorRealMatrix (1, 1) = MatrixCell (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60)
                                      (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60) (-36.9)
generatorRealMatrix (1, 2) = MatrixCell (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60)
                                      (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90) (-84.9)
generatorRealMatrix (1, 3) = MatrixCell (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60)
                                      (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120) (-41.8)
generatorRealMatrix (1, 4) = MatrixCell (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60)
                                      (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150) (-48.5)
generatorRealMatrix (2, 2) = MatrixCell (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90)
                                      (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90) (-68)
generatorRealMatrix (2, 3) = MatrixCell (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90)
                                      (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120) (-92.2)
generatorRealMatrix (2, 4) = MatrixCell (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90)
                                      (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150) (-56.4)
generatorRealMatrix (3, 3) = MatrixCell (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120)
                                      (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120) (-45.7)
generatorRealMatrix (3, 4) = MatrixCell (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120)
                                      (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150) (-81.7)
generatorRealMatrix (4, 4) = MatrixCell (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)
                                      (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150) (-44.5)
generatorRealMatrix (x, y) = MatrixCell (Olig "" 0 0 ) (Olig "" 0 0 ) 0
