{-# LANGUAGE OverloadedStrings #-}

module OligoDesigner.SpecScorer (
    oligoDesignerScoreSpec
) where

import           Bio.Tools.Sequence.OligoDesigner.Scorer (rnaScore, rnaMatrixScore, gcContent,
                                                            oligsGCContentDifference, gcContentScoreBySequence, commonScore)
import           Bio.Tools.Sequence.OligoDesigner.Types  (Olig (..),
                                                          OligLight (..),
                                                          OligSet (..),
                                                          OligSplitting (..), MatrixCell(..))
import           Test.Hspec                              (Spec, describe, it,
                                                          shouldBe)
import Data.Matrix (matrix)

oligoDesignerScoreSpec :: Spec
oligoDesignerScoreSpec = describe "Oligo-Designer score spec" $ do
    scoreRNAFoldingSpec
    scoreRNAFoldingSpec2

    rnaMatrixScoreSpec

    gcContentSpec
    gcContentDifferenceSpec
    gcContentDifferenceForTheSameOligsSpec
    gcContentDifferenceForEmptyOligsSpec

    dnaGCContentScoreSpec
    commonScoreSpec
--       tttt

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
        rnaScore oligs `shouldBe` 13.699997

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
         rnaScore oligs `shouldBe` 1.1000061

rnaMatrixScoreSpec :: Spec
rnaMatrixScoreSpec =
    describe "rnaMatrixScoreSpec" $
    it "" $ do
        let mtx = matrix 4 4 generatorRealMatrix
        let res = rnaMatrixScore mtx
        res `shouldBe` 13.699997
    where
      generatorRealMatrix :: (Int, Int) -> MatrixCell
      generatorRealMatrix (1, 1) = MatrixCell (OligLight "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60))
                                              (OligLight "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60)) (-36.9)
      generatorRealMatrix (1, 2) = MatrixCell (OligLight "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60))
                                              (OligLight "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90)) (-84.9)
      generatorRealMatrix (1, 3) = MatrixCell (OligLight "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60))
                                              (OligLight "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120)) (-41.8)
      generatorRealMatrix (1, 4) = MatrixCell (OligLight "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60))
                                              (OligLight "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)) (-48.5)
      generatorRealMatrix (2, 2) = MatrixCell (OligLight "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90))
                                              (OligLight "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90)) (-68)
      generatorRealMatrix (2, 3) = MatrixCell (OligLight "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90))
                                              (OligLight "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120)) (-92.2)
      generatorRealMatrix (2, 4) = MatrixCell (OligLight "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90))
                                              (OligLight "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)) (-56.4)
      generatorRealMatrix (3, 3) = MatrixCell (OligLight "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120))
                                              (OligLight "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120)) (-45.7)
      generatorRealMatrix (3, 4) = MatrixCell (OligLight "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120))
                                              (OligLight "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)) (-81.7)
      generatorRealMatrix (4, 4) = MatrixCell (OligLight "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150))
                                              (OligLight "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)) (-44.5)
      generatorRealMatrix (_, _) = MatrixCell (OligLight "" (Olig "" 0 0)) (OligLight "" (Olig "" 0 0)) 0

gcContentSpec :: Spec
gcContentSpec = describe "gcContentSpec" $ it "" $ do
    gcContent "GGGGGG" `shouldBe` 100
    gcContent "GGGCCC" `shouldBe` 100
    gcContent "AAAAAA" `shouldBe` 0
    gcContent "TTAAAA" `shouldBe` 0
    gcContent "GGGGCCCCAAAATTTT" `shouldBe` 50
    gcContent "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" `shouldBe` 66.66666666666667

gcContentDifferenceSpec :: Spec
gcContentDifferenceSpec =
    describe "gcContentDifferenceSpec" $
    it "" $ do
        let oligs = OligSet
                        [Olig "TTGATCTTCC" 0 10, Olig "TTATAAGAAA" 10 20] -- 40% & 10%
                        [Olig "TATAAGGAAG" 5 15, Olig "TTGTTTTCT" 15 24]  -- 30% & 22%
                        (OligSplitting [(0, 10), (10, 20)] [(5, 15), (15, 24)])
        let res = oligsGCContentDifference oligs
        res `shouldBe` 30

gcContentDifferenceForTheSameOligsSpec :: Spec
gcContentDifferenceForTheSameOligsSpec =
    describe "gcContentDifferenceForTheSameOligsSpec" $
    it "" $ do
        let oligs = OligSet
                        [Olig "TTGATCTTCC" 0 10, Olig "AAGATGTTCG" 10 20] -- 40% & 10%
                        [Olig "TAGTTCAACC" 5 15, Olig "TTCAACTTGG" 15 24]  -- 30% & 22%
                        (OligSplitting [(0, 10), (10, 20)] [(5, 15), (15, 24)])
        let res = oligsGCContentDifference oligs
        res `shouldBe` 0

gcContentDifferenceForEmptyOligsSpec :: Spec
gcContentDifferenceForEmptyOligsSpec =
    describe "gcContentDifferenceForEmptyOligsSpec" $
    it "" $ do
        let oligs = OligSet [] [] (OligSplitting [(0, 10), (10, 20)] [(5, 15), (15, 24)])
        let res = oligsGCContentDifference oligs
        res `shouldBe` 0

dnaGCContentScoreSpec :: Spec
dnaGCContentScoreSpec =
    describe "dnaGCContentScore" $
    it "" $ do
        let sequ = "GAAGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTACAGCCTGGCAGGTCCCTGTCTCCTGTGCAGCCTCTGGATTCACCTTTAA" ++
                        "CGATTATACCATGCACTGGGTCCAGCTCCAGGGAAGGGCCTGGAGTGGGTCTCAGGTATTAGTTGGAATGGCGGTAG"
        gcContentScoreBySequence "GACATACAAATGACTCAAAGTCCATCTACTCTATCCGCGAGTGTCGGCGACCGCGTAACTATTACGTGCAGGGCTTCACAAAGCATCGGTTCGGCTTTAGCATGGTATCAGCAGAAGCCTGGGAAAGCTCCTAAGTTACTGATCTATAAGGCAAGTGCCCTGGAGAACGGTGTTCCGTCTAGGTTTTCGGGCTCTGGTAGTGGGACCGAGTTCACACTGACAATAAGCAGTCTCCAACCCGATGATTTCGCCACCTACTACTGCCAGCACCTGACCTTCGGACAAGGGACGAGGTTGGAAATCAAA"
            43 `shouldBe` 0.8144094847241223
        gcContentScoreBySequence sequ 64 `shouldBe` 0.9066358024691358
        gcContentScoreBySequence "GGGGGGGGGG" 100 `shouldBe` 1
        gcContentScoreBySequence "AAAAAAAAAA" 100 `shouldBe` 0
        gcContentScoreBySequence "AAAAAGGGGG" 100 `shouldBe` 0.5
        gcContentScoreBySequence "AAAAAGGGGG" 50 `shouldBe` 1
        gcContentScoreBySequence "GGGGGGGGGGGGGGGGGGGG" 50 `shouldBe` 0
        gcContentScoreBySequence "G" 50 `shouldBe` 0
        gcContentScoreBySequence "AAAAAAAAAA" 50 `shouldBe` 0
        gcContentScoreBySequence "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGCGGCACCGCCGCCCTGGG" 71 `shouldBe` 0.993963782696177

commonScoreSpec :: Spec
commonScoreSpec =
    describe "commonScoreSpec" $
    it "" $ do
        let oligs = OligSet
                        [Olig "TTGATCTTCC" 0 10, Olig "TTATAAGAAA" 10 20] -- 40% & 10%
                        [Olig "TATAAGGAAG" 5 15, Olig "TTGTTTTCT" 15 24]  -- 30% & 22%
                        (OligSplitting [(0, 10), (10, 20)] [(5, 15), (15, 24)])
        commonScore 43 oligs `shouldBe` -0.06007751753163892