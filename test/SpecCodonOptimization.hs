{-# LANGUAGE OverloadedStrings #-}

module SpecCodonOptimization where

import           Bio.Protein.AminoAcid                     ()
import           Bio.Tools.Sequence.CodonOptimization.Algo (optimizeAA, optimizeDNA, score,
                                                            scoreCmp)
import           Bio.Tools.Sequence.CodonOptimization.Type (CodonConfig (..), CodonScoreConfig (..))
import           Data.Default                              (def)
import           Data.List                                 (maximumBy,
                                                            minimumBy)
import           Test.Hspec                                (Spec, describe, it,
                                                            shouldBe)

scoreCfg :: CodonScoreConfig
scoreCfg = def

conf :: CodonConfig
conf = def

codonOptimizationSpec :: Spec
codonOptimizationSpec =
    describe "Codon optimization spec" $ do
        optimizeSequence
        optimizeDNASequence
        optimizeShortSequence
        optimizeExtremelyShortSequence
        optimizeSequenceWindow3
        optimizeSequenceInit5
        scoreComparing
        scoreFun
        scoreFunDifferentCodonUsageWeight
        scoreFunDifferentGCWeight
        scoreFunDifferentGCFactor
        scoreFunDifferentGCWindow
        scoreFunDifferentFoldingWeight
        scoreFunDifferentFoldingFactor
        scoreFunDifferentFoldingWindow
        scoreFunWithForbiddenSeq
        scoreFunDifferentForbiddenSeqWeight
        scoreFunDifferentGCDesired

scoreComparing :: Spec
scoreComparing =
    describe "scoreComparing" $
    it "should correct compare by score" $ do
        let optimized = "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGC"
        let vars = ["AGTACTGGT","AGCACTGGT","TCGACAGGT","AGTACAGGT","AGCACAGGT","TCTACTGGC","TCCACTGGC","TCAACTGGC","TCGACTGGC","AGCACCGGC"]
        let cfg = CodonConfig 3 3 (CodonScoreConfig 1 0.5 1.4 40 0.001 2.6 100 1 43)
        let resMin = maximumBy (scoreCmp cfg optimized) vars
        let resMax = minimumBy (scoreCmp cfg optimized) vars
        resMax `shouldBe` "TCGACAGGT"
        resMin `shouldBe` "AGCACCGGC"

optimizeSequence :: Spec
optimizeSequence =
    describe "optimizeSequence" $
    it "should correct optimize sequence" $ do
        let res = optimizeAA conf "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV"
        res `shouldBe` "GCCAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGCGGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGCTGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGCGGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACCTACATCTGCAACGTGAACCACAAGCCTAGCAACACCAAGGTGGACAAGAAGGTG"


optimizeDNASequence :: Spec
optimizeDNASequence =
    describe "optimizeDNASequence" $
    it "should correct optimize amino-acid sequence" $ do
        let res = optimizeDNA conf "GCTAGTACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGCGGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGCTGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGCGGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACCTACATCTGCAACGTGAACCACAAGCCTAGCAACACCAAGGTGGACAAGAAGGTG"
        res `shouldBe` "GCCAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGCGGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGCTGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGCGGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACCTACATCTGCAACGTGAACCACAAGCCTAGCAACACCAAGGTGGACAAGAAGGTG"

optimizeShortSequence :: Spec
optimizeShortSequence =
    describe "optimizeShortSequence" $
    it "should correct optimize short sequence" $ do
        let res = optimizeAA conf "METDTLLLWVLLLWVPGSTG"
        res `shouldBe` "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC"

optimizeExtremelyShortSequence :: Spec
optimizeExtremelyShortSequence =
    describe "optimizeExtremelyShortSequence" $
    it "should correct optimize extremely short sequence" $ do
        let res = optimizeAA conf "METDTLL"
        res `shouldBe` "ATGGAGACCGACACCCTGCTG"

optimizeSequenceWindow3 :: Spec
optimizeSequenceWindow3 =
    describe "optimizeSequenceWindow3" $
    it "should correct optimize sequence with window 3" $ do
        let conf' = CodonConfig 3 3 scoreCfg
        let res = optimizeAA conf' "METDTLLLWVLLLWVPGSTG"
        res `shouldBe` "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCTCTACCGGC"

optimizeSequenceInit5 :: Spec
optimizeSequenceInit5 =
    describe "optimizeSequenceInit5" $
    it "should correct optimize sequence with init param = 5" $ do
        let conf' = CodonConfig 5 1 scoreCfg
        let res = optimizeAA conf' "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV"
        res `shouldBe` "GCCAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGCGGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGCTGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGCGGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACCTACATCTGCAACGTGAACCACAAGCCTAGCAACACCAAGGTGGACAAGAAGGTG"

scoreFun :: Spec
scoreFun =
    describe "scoreFun" $
    it "should correct count score" $ do
        score conf "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` 40.88886033302965
        score conf "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 45.91010855311742
        score conf "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` 66.00968914282373

scoreFunDifferentCodonUsageWeight :: Spec
scoreFunDifferentCodonUsageWeight =
    describe "scoreFunDifferentCodonUsageWeight" $
    it "should correct count score with g_cu=9" $ do
        let conf' = CodonConfig 3 1 (CodonScoreConfig 9 0.5 1.4 40 0.001 2.6 100 1 43)
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` 840.8888603330297
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 748.8141376168824
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` 866.0096891428237

scoreFunDifferentGCWeight :: Spec
scoreFunDifferentGCWeight =
    describe "scoreFunDifferentGCWeight" $
    it "should correct count score with g_gc=0.2" $ do
        let conf' = CodonConfig 3 1 (CodonScoreConfig 1 0.2 1.4 40 0.001 2.6 100 1 43)
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` 75.15554413321186
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 71.08184560102933
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` 85.20387565712949

scoreFunDifferentGCFactor :: Spec
scoreFunDifferentGCFactor =
    describe "scoreFunDifferentGCFactor" $
    it "should correct count score with f_gc=1.7" $ do
        let conf' = CodonConfig 3 1 (CodonScoreConfig 1 0.5 1.7 40 0.001 2.6 100 1 43)
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` -59.640364905662636
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` -20.530452529069947
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` 20.011809536823847

scoreFunDifferentGCWindow :: Spec
scoreFunDifferentGCWindow =
    describe "scoreFunDifferentGCWindow" $
    it "should correct count score with w_gc=10" $ do
        let conf' = CodonConfig 3 1 (CodonScoreConfig 1 0.5 1.4 10 0.001 2.6 100 1 43)
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` 19.575380904184584
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 61.4633338693899
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` 71.60033023641928

scoreFunDifferentFoldingWeight :: Spec
scoreFunDifferentFoldingWeight =
    describe "scoreFunDifferentFoldingWeight" $
    it "should correct count score with weight of folding = 0.03" $ do
        let conf' = CodonConfig 3 1 (CodonScoreConfig 1 0.5 1.4 40 0.03 2.6 100 1 43)
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` -45.11113966697035
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 43.91010855311742
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` -15.990310857176269

scoreFunDifferentFoldingFactor :: Spec
scoreFunDifferentFoldingFactor =
    describe "scoreFunDifferentFoldingFactor" $
    it "should correct count score with folding factor = 4.999" $ do
        let conf' = CodonConfig 3 1 (CodonScoreConfig 1 0.5 1.4 40 0.001 4.999 100 1 43)
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` -4644.11113966697
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 42.91010855311742
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` -4200.990310857176

scoreFunDifferentFoldingWindow :: Spec
scoreFunDifferentFoldingWindow =
    describe "scoreFunDifferentFoldingWindow" $
    it "should correct count score with folding window = 23" $ do
        let conf' = CodonConfig 3 1 (CodonScoreConfig 1 0.5 1.4 40 0.001 2.6 23 1 43)
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` 42.88886033302965
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 45.91010855311742
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` 68.00968914282373

scoreFunWithForbiddenSeq :: Spec
scoreFunWithForbiddenSeq =
    describe "scoreFunWithForbiddenSeq" $
    it "should correct count score with forbidden sequence" $ do
        score conf "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCATACTCCCCCGGC" `shouldBe` -3.4521030544284343
        score conf "GCCAGCGGCGACCGGCGGGTCAAGACCCACACCTGTCCT" `shouldBe` -17.346165749671798
        score conf "ACAGCCAGCGAATAAACCCCGAGGCCGCCGGCGGCCCTAGCGTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGACTACAGCCAAG" `shouldBe` 3.4278532842107197

scoreFunDifferentForbiddenSeqWeight :: Spec
scoreFunDifferentForbiddenSeqWeight =
    describe "scoreFunDifferentForbiddenSeqWeight" $
    it "should correct count score with forbidden sequence weight = 4" $ do
        let conf' = CodonConfig 3 1 (CodonScoreConfig 1 0.5 1.4 40 0.001 2.6 100 4 43)
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCATACTCCCCCGGC" `shouldBe` -153.45210305442845
        score conf' "GCCAGCGGCGACCGGCGGGTCAAGACCCACACCTGTCCT" `shouldBe` -167.3461657496718
        score conf' "ACAGCCAGCGAATAAACCCCGAGGCCGCCGGCGGCCCTAGCGTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGACTACAGCCAAG" `shouldBe` -146.57214671578927
        score conf' "ACAGCCAGCGAATAAACCCCGACGATCGGGCCGCCGGCGAATGACAGCCAAG" `shouldBe` -137.875795016417

scoreFunDifferentGCDesired :: Spec
scoreFunDifferentGCDesired =
    describe "scoreFunDifferentGCDesired" $
    it "should correct count score with gc desired = 60" $ do
        let conf' = CodonConfig 3 1 (CodonScoreConfig 1 0.5 1.4 40 0.001 2.6 100 1 60)
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` 80.83499151084192
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 80.74362252733329
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
         `shouldBe` 96.19662511761598
