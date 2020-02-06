{-# LANGUAGE OverloadedStrings #-}

module OligoDesigner.SpecScorer where

import           Bio.Tools.Sequence.OligoDesigner.Scorer (rnaScore, rnaMatrix, rnaMatrixScore, gcContent,
                                                            gcContentDifference, dnaGCContentScore, commonScore, rebuildMatrix)
import           Bio.Tools.Sequence.OligoDesigner.Types  (Olig (..),
                                                          OligLight (..),
                                                          OligSet (..),
                                                          OligSplitting (..), MatrixCell(..),
                                                          OligsDesignerConfig(..), OligsSplittingConfig(..))
import           Test.Hspec                              (Spec, describe, it,
                                                          shouldBe)
import Data.Matrix (matrix)
import Bio.Tools.Sequence.OligoDesigner.Utils (assemble)
import Debug.Trace (trace)
import Data.Text (toUpper)
import Data.Default (def)
import Bio.Tools.Sequence.CodonOptimization (CodonOptimizationConfig(..))
import Bio.Tools.Sequence.CodonOptimization.Types (Organism(..), defaultForbiddenRegexp)
import Bio.NucleicAcid.Nucleotide (DNA)

oligoDesignerScoreSpec :: Spec
oligoDesignerScoreSpec = describe "Oligo-Designer score spec" $ do
    scoreRNAFoldingSpec
    scoreRNAFoldingSpec2

    rnaMatrixSpec
    rnaMatrixSimpleSpec
    rnaMatrixScoreSpec

    gcContentSpec
    gcContentDifferenceSpec
    gcContentDifferenceForTheSameOligsSpec
    gcContentDifferenceForEmptyOligsSpec

    dnaGCContentScoreSpec
    commonScoreSpec

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
        generator (1, 1) = MatrixCell (OligLight (dnaToStr "AATGGC") (Olig "AATGGC" 0 6))
                                      (OligLight (dnaToStr "AATGGC") (Olig "AATGGC" 0 6)) 0
        generator (1, 2) = MatrixCell (OligLight (dnaToStr "AATGGC") (Olig "AATGGC" 0 6))
                                      (OligLight (dnaToStr "CCGCTG") (Olig "CCGCTG" 3 9)) (-2)
        generator (1, 3) = MatrixCell (OligLight (dnaToStr "AATGGC") (Olig "AATGGC" 0 6))
                                      (OligLight (dnaToStr "GACTGA") (Olig "GACTGA" 6 12)) 0
        generator (1, 4) = MatrixCell (OligLight (dnaToStr "AATGGC") (Olig "AATGGC" 0 6))
                                      (OligLight (dnaToStr "ACTTTC") (Olig "ACTTTC" 9 15)) 0
        generator (2, 2) = MatrixCell (OligLight (dnaToStr "CCGCTG") (Olig "CCGCTG" 3 9))
                                      (OligLight (dnaToStr "CCGCTG") (Olig "CCGCTG" 3 9)) (-2.1)
        generator (2, 3) = MatrixCell (OligLight (dnaToStr "CCGCTG") (Olig "CCGCTG" 3 9))
                                      (OligLight (dnaToStr "GACTGA") (Olig "GACTGA" 6 12)) 0
        generator (2, 4) = MatrixCell (OligLight (dnaToStr "CCGCTG") (Olig "CCGCTG" 3 9))
                                      (OligLight (dnaToStr "ACTTTC") (Olig "ACTTTC" 9 15)) 0
        generator (3, 3) = MatrixCell (OligLight (dnaToStr "GACTGA") (Olig "GACTGA" 6 12))
                                      (OligLight (dnaToStr "GACTGA") (Olig "GACTGA" 6 12)) 0
        generator (3, 4) = MatrixCell (OligLight (dnaToStr "GACTGA") (Olig "GACTGA" 6 12))
                                      (OligLight (dnaToStr "ACTTTC") (Olig "ACTTTC" 9 15)) 0
        generator (4, 4) = MatrixCell (OligLight (dnaToStr "ACTTTC") (Olig "ACTTTC" 9 15))
                                      (OligLight (dnaToStr "ACTTTC") (Olig "ACTTTC" 9 15)) 0
        generator (x, y) = MatrixCell (OligLight "" (Olig "" 0 0 )) (OligLight "" (Olig "" 0 0)) 0

rnaMatrixScoreSpec :: Spec
rnaMatrixScoreSpec =
    describe "rnaMatrixScoreSpec" $
    it "" $ do
        let mtx = matrix 4 4 generatorRealMatrix
        let res = rnaMatrixScore mtx
        res `shouldBe` -55.299995

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
        let res = gcContentDifference oligs
        res `shouldBe` 30

gcContentDifferenceForTheSameOligsSpec :: Spec
gcContentDifferenceForTheSameOligsSpec =
    describe "gcContentDifferenceForTheSameOligsSpec" $
    it "" $ do
        let oligs = OligSet
                        [Olig "TTGATCTTCC" 0 10, Olig "AAGATGTTCG" 10 20] -- 40% & 10%
                        [Olig "TAGTTCAACC" 5 15, Olig "TTCAACTTGG" 15 24]  -- 30% & 22%
                        (OligSplitting [(0, 10), (10, 20)] [(5, 15), (15, 24)])
        let res = gcContentDifference oligs
        res `shouldBe` 0

gcContentDifferenceForEmptyOligsSpec :: Spec
gcContentDifferenceForEmptyOligsSpec =
    describe "gcContentDifferenceForEmptyOligsSpec" $
    it "" $ do
        let oligs = OligSet [] [] (OligSplitting [(0, 10), (10, 20)] [(5, 15), (15, 24)])
        let res = gcContentDifference oligs
        res `shouldBe` 0

dnaGCContentScoreSpec :: Spec
dnaGCContentScoreSpec =
    describe "dnaGCContentScore" $
    it "" $ do
        let sequ = "GAAGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTACAGCCTGGCAGGTCCCTGTCTCCTGTGCAGCCTCTGGATTCACCTTTAA" ++
                        "CGATTATACCATGCACTGGGTCCAGCTCCAGGGAAGGGCCTGGAGTGGGTCTCAGGTATTAGTTGGAATGGCGGTAG"
        dnaGCContentScore sequ 64 `shouldBe` 0.9066358024691358
        dnaGCContentScore "GGGGGGGGGG" 100 `shouldBe` 1
        dnaGCContentScore "AAAAAAAAAA" 100 `shouldBe` 0
        dnaGCContentScore "AAAAAGGGGG" 100 `shouldBe` 0.5
        dnaGCContentScore "AAAAAGGGGG" 50 `shouldBe` 1
        dnaGCContentScore "GGGGGGGGGGGGGGGGGGGG" 50 `shouldBe` 0
        dnaGCContentScore "G" 50 `shouldBe` 0
        dnaGCContentScore "AAAAAAAAAA" 50 `shouldBe` 0
        dnaGCContentScore "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGCGGCACCGCCGCCCTGGG" 71 `shouldBe` 0.993963782696177 --71%

commonScoreSpec :: Spec
commonScoreSpec =
    describe "commonScoreSpec" $
    it "" $ do
        let oligs = OligSet
                        [Olig "TTGATCTTCC" 0 10, Olig "TTATAAGAAA" 10 20] -- 40% & 10%
                        [Olig "TATAAGGAAG" 5 15, Olig "TTGTTTTCT" 15 24]  -- 30% & 22%
                        (OligSplitting [(0, 10), (10, 20)] [(5, 15), (15, 24)])
        let codonConf = CodonOptimizationConfig CHO 3 1 1 0.5 1.4 40 0.001 2.6 100 1 43 defaultForbiddenRegexp
        commonScore (OligsDesignerConfig codonConf def 1 0 0 0 1) oligs `shouldBe` -4.900000095367432
        commonScore (OligsDesignerConfig codonConf def 0 1 0 0 1) oligs `shouldBe` 30
        commonScore (OligsDesignerConfig codonConf def 0 0 1 0 1) oligs `shouldBe` 0.5813953488372092
        commonScore (OligsDesignerConfig codonConf def 1 1 0 0 1) oligs `shouldBe` 25.09999990463257
        commonScore (OligsDesignerConfig codonConf def 1 1 1 0 1) oligs `shouldBe` 25.681395253469777
        commonScore (OligsDesignerConfig codonConf def 0 1 1 0 1) oligs `shouldBe` 30.5813953488372092
        commonScore (OligsDesignerConfig codonConf def 0 0 0 0 1) oligs `shouldBe` 0

generatorRealMatrix :: (Int, Int) -> MatrixCell
generatorRealMatrix (1, 1) = MatrixCell (OligLight (dnaToStr "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC") (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60))
                                        (OligLight (dnaToStr "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC") (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60)) (-36.9)
generatorRealMatrix (1, 2) = MatrixCell (OligLight (dnaToStr "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC") (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60))
                                        (OligLight (dnaToStr "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG") (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90)) (-84.9)
generatorRealMatrix (1, 3) = MatrixCell (OligLight (dnaToStr "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC") (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60))
                                        (OligLight (dnaToStr "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC") (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120)) (-41.8)
generatorRealMatrix (1, 4) = MatrixCell (OligLight (dnaToStr "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC") (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60))
                                        (OligLight (dnaToStr "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC") (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)) (-48.5)
generatorRealMatrix (2, 2) = MatrixCell (OligLight (dnaToStr "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG") (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90))
                                        (OligLight (dnaToStr "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG") (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90)) (-68)
generatorRealMatrix (2, 3) = MatrixCell (OligLight (dnaToStr "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG") (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90))
                                        (OligLight (dnaToStr "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC") (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120)) (-92.2)
generatorRealMatrix (2, 4) = MatrixCell (OligLight (dnaToStr "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG") (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90))
                                        (OligLight (dnaToStr "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC") (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)) (-56.4)
generatorRealMatrix (3, 3) = MatrixCell (OligLight (dnaToStr "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC") (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120))
                                        (OligLight (dnaToStr "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC") (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120)) (-45.7)
generatorRealMatrix (3, 4) = MatrixCell (OligLight (dnaToStr "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC") (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120))
                                        (OligLight (dnaToStr "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC") (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)) (-81.7)
generatorRealMatrix (4, 4) = MatrixCell (OligLight (dnaToStr "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC") (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150))
                                        (OligLight (dnaToStr "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC") (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)) (-44.5)
generatorRealMatrix (x, y) = MatrixCell (OligLight "" (Olig "" 0 0)) (OligLight "" (Olig "" 0 0)) 0


rebuildMatrixSpec :: Spec
rebuildMatrixSpec =
    describe "rebuildMatrixSpec" $
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
        let mtx = rnaMatrix oligs
        let newOligs =
                OligSet
                    [ Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGA" 0 60
                    , Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120
                    ]
                    [ Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTTCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90
                    , Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150
                    ]
                    coords
        let newMtx = rebuildMatrix mtx newOligs
        newMtx `shouldBe` matrix 4 4 generator
  where
    generator :: (Int, Int) -> MatrixCell
    generator (1, 1) = MatrixCell (OligLight (dnaToStr "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGA") (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGA" 0 60))
                                  (OligLight (dnaToStr "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGA") (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGA" 0 60)) (-40.1)
    generator (1, 2) = MatrixCell (OligLight (dnaToStr "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGA") (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGA" 0 60))
                                  (OligLight (dnaToStr "CTTCACCAGGCAGCCCAGGGCGGCGGTTCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG") (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTTCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90)) (-79.3)
    generator (1, 3) = MatrixCell (OligLight (dnaToStr "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGA") (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGA" 0 60))
                                  (OligLight (dnaToStr "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC") (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120)) (-40.1)
    generator (1, 4) = MatrixCell (OligLight (dnaToStr "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGA") (Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGA" 0 60))
                                  (OligLight (dnaToStr "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC") (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)) (-46.1)
    generator (2, 2) = MatrixCell (OligLight (dnaToStr "CTTCACCAGGCAGCCCAGGGCGGCGGTTCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG") (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTTCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90))
                                  (OligLight (dnaToStr "CTTCACCAGGCAGCCCAGGGCGGCGGTTCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG") (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTTCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90)) (-60.7)
    generator (2, 3) = MatrixCell (OligLight (dnaToStr "CTTCACCAGGCAGCCCAGGGCGGCGGTTCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG") (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTTCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90))
                                  (OligLight (dnaToStr "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC") (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120)) (-85.5)
    generator (2, 4) = MatrixCell (OligLight (dnaToStr "CTTCACCAGGCAGCCCAGGGCGGCGGTTCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG") (Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTTCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90))
                                  (OligLight (dnaToStr "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC") (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)) (-54.0)
    generator (3, 3) = MatrixCell (OligLight (dnaToStr "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC") (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120))
                                  (OligLight (dnaToStr "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC") (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120)) (-45.7)
    generator (3, 4) = MatrixCell (OligLight (dnaToStr "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC") (Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120))
                                  (OligLight (dnaToStr "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC") (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)) (-81.7)
    generator (4, 4) = MatrixCell (OligLight (dnaToStr "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC") (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150))
                                  (OligLight (dnaToStr "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC") (Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150)) (-44.5)
    generator (x, y) = MatrixCell (OligLight "" (Olig "" 0 0)) (OligLight "" (Olig "" 0 0)) 0
    
dnaToStr :: [DNA] -> String
dnaToStr = show