{-# LANGUAGE OverloadedStrings #-}

module OligoDesigner.SpecOligoDesignerScorer where

import           Bio.Tools.Sequence.OligoDesigner.Scorer (rnaScore, rnaMatrix)
import           Bio.Tools.Sequence.OligoDesigner.Types  (Olig (..),
                                                          OligSet (..),
                                                          OligSplitting (..), MatrixCell(..))
import           Test.Hspec                              (Spec, describe, it,
                                                          shouldBe)
import Data.Matrix (matrix)
import Bio.Tools.Sequence.OligoDesigner.Utils (assemble)

oligoDesignerScoreSpec :: Spec
oligoDesignerScoreSpec = describe "Oligo-Designer score spec" $ do
--    scoreRNAFoldingSpec
--    scoreRNAFoldingSpec2
    rnaMatrixSpec
    rnaMatrixSimpleSpec

scoreRNAFoldingSpec :: Spec
scoreRNAFoldingSpec =
    describe "scoreRNAFoldingSpec" $
    it "should correct score oligs" $ do
        let fwd =
                [ Olig "TCTCTGGCCTAACTGGCCGGTACCTGAGCTCAGGGTTCCTGGAAGGGAATTATT" 0 54
                , Olig "CTGAGAAAAGGAGCTTAGCTATGACTCATCCGGAAGCTCAGATTCCGGGAATCC" 54 108
                , Olig "CTACTTTCTTAGAAAATTGTATTAGTCATCGCTATTACCATGATCCTTCTGGGA" 108 162
                , Olig "ATTCTGGTGTTCCTAGAAGAGGGATAGCGGTTTGACTCACGGGGATGATATCAA" 162 216
                ]
        let rvsd =
                [ Olig "GAGTCATAGCTAAGCTCCTTTTCTCAGAATAATTCCCTTCCAGGAACCCTGAGC" 27 81
                , Olig "GACTAATACAATTTTCTAAGAAAGTAGGGATTCCCGGAATCTGAGCTTCCGGAT" 81 135
                , Olig "CTATCCCTCTTCTAGGAACACCAGAATTCCCAGAAGGATCATGGTAATAGCGAT" 135 189
                , Olig "ATATACCCTCTAGAGTCTAGATCTTGATATCATCCCCGTGAGTCAAACCG" 189 239
                ]
        let oligSet = OligSet fwd rvsd (OligSplitting [] [])
        rnaScore oligSet `shouldBe` -7.6000023

scoreRNAFoldingSpec2 :: Spec
scoreRNAFoldingSpec2 =
    describe "scoreRNAFoldingSpec2" $
    it "should correct score oligs" $ do
        let fwd =
                [ Olig "GAAGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTACAGCCTGGCAGGTCCCTG" 0 0
                , Olig "TCTCCTGTGCAGCCTCTGGATTCACCTTTAACGATTATACCATGCACTGGGTCC" 0 0
                , Olig "AGCTCCAGGGAAGGGCCTGGAGTGGGTCTCAGGTATTAGTTGGAATGGCGGTAG" 0 0
                , Olig "GGCTATGCGGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAG" 0 0
                , Olig "CCCTGTATCTGCAAATGAACAGTCTGAGAGCTGAGGACACGGCCTTGTATTACT" 0 0
                , Olig "AAAAGATAATAGTGGCTACGGTCACTACTACTACGGAATGGACATCTGGGGCCA" 0 0
                ]
        let rvsd =
                [ Olig "GTGAATCCAGAGGCTGCACAGGAGAGTCTCAGGGACCTGCCAGGCTGTACCACG" 0 0
                , Olig "CCACTCCAGGCCCTTCCCTGGAGCTTGCCGGACCCAGTGCATGGTATAATCGTT" 0 0
                , Olig "GGCCCTTCACAGAGTCCGCATAGCCTATACTACCGCCATTCCAACTAATACCTG" 0 0
                , Olig "AGACTGTTCATTTGCAGATACAGGGAGTTCTTGGCGTTGTCTCTGGAGATGGTG" 0 0
                , Olig "GTGACCGTAGCCACTATTATCTTTTGCACAGTAATACAAGGCCGTGTCCTCAGC" 0 0
                , Olig "TGAACTGACGGTGACCGTGGTCCCTTGGCCCCAGATGTCCATTCCGTAGT" 0 0
                ]
        let oligSet = OligSet fwd rvsd (OligSplitting [] [])
        rnaScore oligSet `shouldBe` 5.5999985

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
        res `shouldBe` matrix 4 4 generator
      where
        generator :: (Int, Int) -> MatrixCell
        generator (1, 1) = MatrixCell (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60)
                                      (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60) (-36.9)
        generator (1, 2) = MatrixCell (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60)
                                      (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90) (-84.9)
        generator (1, 3) = MatrixCell (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60)
                                      (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120) (-41.8)
        generator (1, 4) = MatrixCell (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60)
                                      (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150) (-48.5)
        generator (2, 2) = MatrixCell (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90)
                                      (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90) (-68)
        generator (2, 3) = MatrixCell (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90)
                                      (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120) (-92.2)
        generator (2, 4) = MatrixCell (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90)
                                      (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150) (-56.4)
        generator (3, 3) = MatrixCell (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120)
                                      (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120) (-45.7)
        generator (3, 4) = MatrixCell (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120)
                                      (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150) (-81.7)
        generator (4, 4) = MatrixCell (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)
                                      (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150) (-44.5)
        generator (x, y) = MatrixCell (Olig "" 0 0 ) (Olig "" 0 0 ) 0

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

