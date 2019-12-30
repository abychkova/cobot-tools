{-# LANGUAGE OverloadedStrings #-}

module SpecOligoDesignerScorer where

import Bio.Tools.Sequence.OligoDesigner.Scorer (rnaCofoldScore, gcScore)
import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligSet(..))
import Test.Hspec (Expectation, Spec, describe, it, shouldBe, shouldSatisfy)
import Test.Hspec.QuickCheck (modifyMaxSize, prop)
import Test.QuickCheck (Gen, Property, elements, forAll)

oligoDesignerScoreSpec :: Spec
oligoDesignerScoreSpec = describe "Oligo-Designer score spec" $ do
    scoreRNAFoldingSpec
    scoreRNAFoldingSpec2
    gcScoreSpec

scoreRNAFoldingSpec :: Spec
scoreRNAFoldingSpec =
    describe "scoreRNAFoldingSpec" $
    it "should correct score oligs" $ do
        let forward =
                [ Olig "TCTCTGGCCTAACTGGCCGGTACCTGAGCTCAGGGTTCCTGGAAGGGAATTATT" 0 54
                , Olig "CTGAGAAAAGGAGCTTAGCTATGACTCATCCGGAAGCTCAGATTCCGGGAATCC" 54 108
                , Olig "CTACTTTCTTAGAAAATTGTATTAGTCATCGCTATTACCATGATCCTTCTGGGA" 108 162
                , Olig "ATTCTGGTGTTCCTAGAAGAGGGATAGCGGTTTGACTCACGGGGATGATATCAA" 162 216
                ]
        let reversed =
                [ Olig "GAGTCATAGCTAAGCTCCTTTTCTCAGAATAATTCCCTTCCAGGAACCCTGAGC" 27 81
                , Olig "GACTAATACAATTTTCTAAGAAAGTAGGGATTCCCGGAATCTGAGCTTCCGGAT" 81 135
                , Olig "CTATCCCTCTTCTAGGAACACCAGAATTCCCAGAAGGATCATGGTAATAGCGAT" 135 189
                , Olig "ATATACCCTCTAGAGTCTAGATCTTGATATCATCCCCGTGAGTCAAACCG" 189 239
                ]
        let oligSet = OligSet forward reversed
        rnaCofoldScore oligSet `shouldBe` -7.6000023

scoreRNAFoldingSpec2 :: Spec
scoreRNAFoldingSpec2 =
    describe "scoreRNAFoldingSpec2" $
    it "should correct score oligs" $ do
        let forward =
                [ Olig "GAAGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTACAGCCTGGCAGGTCCCTG" 0 0
                , Olig "TCTCCTGTGCAGCCTCTGGATTCACCTTTAACGATTATACCATGCACTGGGTCC" 0 0
                , Olig "AGCTCCAGGGAAGGGCCTGGAGTGGGTCTCAGGTATTAGTTGGAATGGCGGTAG" 0 0
                , Olig "GGCTATGCGGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAG" 0 0
                , Olig "CCCTGTATCTGCAAATGAACAGTCTGAGAGCTGAGGACACGGCCTTGTATTACT" 0 0
                , Olig "AAAAGATAATAGTGGCTACGGTCACTACTACTACGGAATGGACATCTGGGGCCA" 0 0
                ]
        let reversed =
                [ Olig "GTGAATCCAGAGGCTGCACAGGAGAGTCTCAGGGACCTGCCAGGCTGTACCACG" 0 0
                , Olig "CCACTCCAGGCCCTTCCCTGGAGCTTGCCGGACCCAGTGCATGGTATAATCGTT" 0 0
                , Olig "GGCCCTTCACAGAGTCCGCATAGCCTATACTACCGCCATTCCAACTAATACCTG" 0 0
                , Olig "AGACTGTTCATTTGCAGATACAGGGAGTTCTTGGCGTTGTCTCTGGAGATGGTG" 0 0
                , Olig "GTGACCGTAGCCACTATTATCTTTTGCACAGTAATACAAGGCCGTGTCCTCAGC" 0 0
                , Olig "TGAACTGACGGTGACCGTGGTCCCTTGGCCCCAGATGTCCATTCCGTAGT" 0 0
                ]
        let oligSet = OligSet forward reversed
        rnaCofoldScore oligSet `shouldBe` 5.5999985
        
gcScoreSpec :: Spec
gcScoreSpec = 
   describe "" $ it "" $ do
       let sequ = "GAAGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTACAGCCTGGCAGGTCCCTGTCTCCTGTGCAGCCTCTGGATTCACCTTTAA" ++ 
                  "CGATTATACCATGCACTGGGTCCAGCTCCAGGGAAGGGCCTGGAGTGGGTCTCAGGTATTAGTTGGAATGGCGGTAG"
       let target = 64
       gcScore sequ target `shouldBe` 0.009066358024691357