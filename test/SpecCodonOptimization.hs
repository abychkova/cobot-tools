{-# LANGUAGE OverloadedStrings #-}

module SpecCodonOptimization where

import           Bio.NucleicAcid.Nucleotide.Type                (DNA)
import           Bio.Protein.AminoAcid                          ()
import           Bio.Protein.AminoAcid.Type                     (AA)
import           Bio.Tools.Sequence.CodonOptimization.Algo      (optimizeAA,
                                                                 optimizeDNA,
                                                                 score,
                                                                 scoreCmp,
                                                                 scoreSequence)
import           Bio.Tools.Sequence.CodonOptimization.Constants (ak2Codon)
import           Bio.Tools.Sequence.CodonOptimization.Types     (CodonScoreConfig(..),
                                                                 Organism (..),
                                                                 defaultForbiddenRegexp)
import           Data.List                                      (foldl',
                                                                 maximumBy,
                                                                 minimumBy)
import           Data.Map                                       as Map (lookup)
import           Data.Maybe                                     (fromMaybe)
import           System.Random
import           Test.Hspec                                     (Expectation,
                                                                 Spec, describe,
                                                                 it, shouldBe,
                                                                 shouldSatisfy)

confHuman :: CodonScoreConfig
confHuman = CodonScoreConfig Human 3 1 1 0.5 1.4 40 0.001 2.6 100 1 60 defaultForbiddenRegexp

confEColi :: CodonScoreConfig
confEColi = CodonScoreConfig EColi 3 1 1 0.5 1.4 40 0.001 2.6 100 1 50 defaultForbiddenRegexp

confCHO :: CodonScoreConfig
confCHO = CodonScoreConfig CHO 3 1 1 0.5 1.4 40 0.001 2.6 100 1 43 defaultForbiddenRegexp

scoreHuman :: [DNA] -> Double
scoreHuman = scoreSequence confHuman

scoreCHO :: [DNA] -> Double
scoreCHO = scoreSequence confCHO

scoreEColi :: [DNA] -> Double
scoreEColi = scoreSequence confEColi

assertScoreBecomeBetter :: CodonScoreConfig -> [DNA] -> [DNA] -> Expectation
assertScoreBecomeBetter conf res initial = scoreSequence conf initial `shouldSatisfy` (< scoreSequence conf res)

assertScoreByHumanBetterThan :: [DNA] -> [DNA] -> Expectation
assertScoreByHumanBetterThan res initial = scoreHuman initial `shouldSatisfy` (< scoreHuman res)

assertScoreByCHOBetterThan :: [DNA] -> [DNA] -> Expectation
assertScoreByCHOBetterThan res initial = scoreCHO initial `shouldSatisfy` (< scoreCHO res)

assertScoreByEColiBetterThan :: [DNA] -> [DNA] -> Expectation
assertScoreByEColiBetterThan res initial = scoreEColi initial `shouldSatisfy` (< scoreEColi res)

toRandomNKSequ :: [AA] -> IO [DNA]
toRandomNKSequ initial = foldl' (++) [] <$> mapM toRandomCodon initial

toRandomCodon :: AA -> IO [DNA]
toRandomCodon ak = do
    let codons = fromMaybe [] (Map.lookup ak ak2Codon)
    i <- randomRIO (0, length codons - 1)
    return $ codons !! i

codonOptimizationSpec :: Spec
codonOptimizationSpec =
    describe "Codon optimization spec" $ do
        optimizeSequence
        optimizeSequenceForEColi
        optimizeSequenceForCHO
        optimizeDNASequence
        optimizeShortSequence
        optimizeExtremelyShortSequence
        optimizeSequenceWindow3
        optimizeSequenceInit5
        scoreComparing
        scoreFun
        scoreFunEColi
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
        let cfg = CodonScoreConfig Human 3 3 1 0.5 1.4 40 0.001 2.6 100 1 43 defaultForbiddenRegexp
        let resMin = maximumBy (scoreCmp cfg optimized) vars
        let resMax = minimumBy (scoreCmp cfg optimized) vars
        resMax `shouldBe` "TCGACAGGT"
        resMin `shouldBe` "AGCACCGGC"

optimizeSequence :: Spec
optimizeSequence =
    describe "optimizeSequence" $
    it "should correct optimize sequence" $ do
        let ak = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV"
        nk <- toRandomNKSequ ak
        let res = optimizeAA confHuman ak

        res `shouldBe` "GCCAGCACCAAGGGCCCCAGCGTGTTCCCTCTGGCCCCCAGCAGCAAGAGCACCAGCGGCGGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCCGAGCCCGTGACCGTGAGCTGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGCGGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACCTACATCTGCAACGTGAACCACAAGCCCAGCAACACCAAGGTGGACAAGAAGGTG"
        res `assertScoreByHumanBetterThan` nk

optimizeSequenceForCHO :: Spec
optimizeSequenceForCHO =
    describe "optimizeSequenceForCHO" $
    it "should correct optimize sequence for CHO" $ do
        let ak = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV"
        nk <- toRandomNKSequ ak
        let res = optimizeAA confCHO ak

        res `shouldBe` "GCCAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTTCTTCTAAGTCTACCTCTGGCGGCACCGCCGCCCTGGGCTGTCTGGTGAAGGATTACTTTCCTGAGCCTGTGACCGTGTCTTGGAACTCTGGCGCCCTGACCTCTGGCGTGCACACCTTTCCTGCCGTGCTGCAGTCTTCTGGCCTGTACTCTCTGTCTTCTGTGGTGACCGTGCCTTCTTCTTCTCTGGGCACCCAGACCTACATCTGCAACGTGAACCACAAGCCTTCTAACACCAAGGTGGACAAGAAGGTG"
        res `assertScoreByCHOBetterThan` nk

optimizeSequenceForEColi :: Spec
optimizeSequenceForEColi =
    describe "optimizeSequenceForEColi" $
    it "should correct optimize sequence for Ecoli" $ do
        let ak = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV"
        nk <- toRandomNKSequ ak
        let res = optimizeAA confEColi ak

        res `shouldBe` "GCGAGCACCAAAGGCCCGAGCGTGTTTCCGCTGGCGCCGAGCAGCAAAAGCACCAGCGGCGGCACCGCGGCGCTGGGCTGCCTGGTGAAAGATTATTTTCCGGAACCGGTGACCGTGAGCTGGAACAGCGGCGCGCTGACCAGCGGCGTGCATACCTTTCCGGCGGTGCTGCAGAGCAGCGGCCTGTATAGCCTGAGCAGCGTGGTGACCGTGCCGAGCAGCAGCCTGGGCACCCAGACCTATATTTGCAACGTGAACCATAAACCGAGCAACACCAAAGTGGATAAAAAAGTG"
        res `assertScoreByEColiBetterThan` nk

optimizeDNASequence :: Spec
optimizeDNASequence =
    describe "optimizeDNASequence" $
    it "should correct optimize amino-acid sequence" $ do
        let nk = "GCTAGCACGAAAGGCCCTTCAGTATTCCCCCTCGCACCGTCGAGTAAGTCCACGTCGGGTGGGACGGCGGCTCTAGGATGCTTAGTTAAGGACTATTTTCCAGAGCCTGTCACAGTGTCGTGGAACAGTGGTGCTTTAACCAGCGGTGTCCACACCTTTCCTGCCGTTTTACAAAGTAGTGGTCTTTATTCCCTATCGAGCGTCGTTACGGTTCCCAGTTCGAGTTTGGGGACACAGACATACATTTGTAACGTAAACCACAAACCCTCTAACACGAAAGTCGATAAGAAAGTC"
        let res = optimizeDNA confHuman nk

        res `shouldBe` "GCCAGCACCAAGGGCCCCAGCGTGTTCCCTCTGGCCCCCAGCAGCAAGAGCACCAGCGGCGGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCCGAGCCCGTGACCGTGAGCTGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGCGGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACCTACATCTGCAACGTGAACCACAAGCCCAGCAACACCAAGGTGGACAAGAAGGTG"
        res `assertScoreByHumanBetterThan` nk

optimizeShortSequence :: Spec
optimizeShortSequence =
    describe "optimizeShortSequence" $
    it "should correct optimize short sequence" $ do
        let ak = "METDTLLLWVLLLWVPGSTG"
        nk <- toRandomNKSequ ak
        let res = optimizeAA confHuman ak

        res `shouldBe` "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCCGGCAGCACCGGC"
        res `assertScoreByHumanBetterThan` nk

optimizeExtremelyShortSequence :: Spec
optimizeExtremelyShortSequence =
    describe "optimizeExtremelyShortSequence" $
    it "should correct optimize extremely short sequence" $ do
        let ak = "METDTLL"
        nk <- toRandomNKSequ ak
        let res = optimizeAA confHuman ak

        res `shouldBe` "ATGGAGACCGACACCCTGCTG"
        res `assertScoreByHumanBetterThan` nk

optimizeSequenceWindow3 :: Spec
optimizeSequenceWindow3 =
    describe "optimizeSequenceWindow3" $
    it "should correct optimize sequence with window 3" $ do
        let conf' = CodonScoreConfig Human 3 3 1 0.5 1.4 40 0.001 2.6 100 1 60 defaultForbiddenRegexp
        let ak = "METDTLLLWVLLLWVPGSTG"
        nk <- toRandomNKSequ ak
        let res = optimizeAA conf' ak

        res `shouldBe` "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC"
        assertScoreBecomeBetter conf' res nk

optimizeSequenceInit5 :: Spec
optimizeSequenceInit5 =
    describe "optimizeSequenceInit5" $
    it "should correct optimize sequence with init param = 5" $ do
        let conf' = CodonScoreConfig Human 5 1 1 0.5 1.4 40 0.001 2.6 100 1 60 defaultForbiddenRegexp
        let ak = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV"
        nk <- toRandomNKSequ ak
        let res = optimizeAA conf' ak

        res `shouldBe` "GCCAGCACCAAGGGCCCCAGCGTGTTCCCTCTGGCCCCCAGCAGCAAGAGCACCAGCGGCGGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCCGAGCCCGTGACCGTGAGCTGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGCGGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACCTACATCTGCAACGTGAACCACAAGCCCAGCAACACCAAGGTGGACAAGAAGGTG"
        assertScoreBecomeBetter conf' res nk

scoreFun :: Spec
scoreFun =
    describe "scoreFun" $
    it "should correct count score" $ do
        score confHuman "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` 80.83499151084192
        score confHuman "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 80.74362252733329
        score confHuman "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` 96.19662511761598

scoreFunEColi :: Spec
scoreFunEColi =
    describe "scoreFunEcoli" $
    it "should correct count score for Ecoli" $ do
        score confEColi "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` 58.913612822848826
        score confEColi "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 23.025024455733913
        score confEColi "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` 29.501417850070677

scoreFunDifferentCodonUsageWeight :: Spec
scoreFunDifferentCodonUsageWeight =
    describe "scoreFunDifferentCodonUsageWeight" $
    it "should correct count score with g_cu=9" $ do
        let conf' = CodonScoreConfig Human 3 1 9 0.5 1.4 40 0.001 2.6 100 1 43 defaultForbiddenRegexp
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` 840.8888603330297
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 748.8141376168824
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` 866.0096891428237

scoreFunDifferentGCWeight :: Spec
scoreFunDifferentGCWeight =
    describe "scoreFunDifferentGCWeight" $
    it "should correct count score with g_gc=0.2" $ do
        let conf' = CodonScoreConfig Human 3 1 1 0.2 1.4 40 0.001 2.6 100 1 43 defaultForbiddenRegexp
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` 75.15554413321186
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 71.08184560102933
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` 85.20387565712949

scoreFunDifferentGCFactor :: Spec
scoreFunDifferentGCFactor =
    describe "scoreFunDifferentGCFactor" $
    it "should correct count score with f_gc=1.7" $ do
        let conf' = CodonScoreConfig Human 3 1 1 0.5 1.7 40 0.001 2.6 100 1 43 defaultForbiddenRegexp
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` -59.640364905662636
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` -20.530452529069947
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` 20.011809536823847

scoreFunDifferentGCWindow :: Spec
scoreFunDifferentGCWindow =
    describe "scoreFunDifferentGCWindow" $
    it "should correct count score with w_gc=10" $ do
        let conf' = CodonScoreConfig Human 3 1 1 0.5 1.4 10 0.001 2.6 100 1 43 defaultForbiddenRegexp
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` 19.575380904184584
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 61.4633338693899
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` 71.60033023641928

scoreFunDifferentFoldingWeight :: Spec
scoreFunDifferentFoldingWeight =
    describe "scoreFunDifferentFoldingWeight" $
    it "should correct count score with weight of folding = 0.03" $ do
        let conf' = CodonScoreConfig Human 3 1 1 0.5 1.4 40 0.03 2.6 100 1 43 defaultForbiddenRegexp
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` -45.11113966697035
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 43.91010855311742
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` -15.990310857176269

scoreFunDifferentFoldingFactor :: Spec
scoreFunDifferentFoldingFactor =
    describe "scoreFunDifferentFoldingFactor" $
    it "should correct count score with folding factor = 4.999" $ do
        let conf' = CodonScoreConfig Human 3 1 1 0.5 1.4 40 0.001 4.999 100 1 43 defaultForbiddenRegexp
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` -4644.11113966697
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 42.91010855311742
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` -4200.990310857176

scoreFunDifferentFoldingWindow :: Spec
scoreFunDifferentFoldingWindow =
    describe "scoreFunDifferentFoldingWindow" $
    it "should correct count score with folding window = 23" $ do
        let conf' = CodonScoreConfig Human 3 1 1 0.5 1.4 40 0.001 2.6 23 1 43 defaultForbiddenRegexp
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` 42.88886033302965
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 45.91010855311742
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
            `shouldBe` 68.00968914282373

scoreFunWithForbiddenSeq :: Spec
scoreFunWithForbiddenSeq =
    describe "scoreFunWithForbiddenSeq" $
    it "should correct count score with forbidden sequence" $ do
        let confWithoutForbidden = CodonScoreConfig Human 3 1 1 0.5 1.4 40 0.001 2.6 23 1 43 []
        score confWithoutForbidden "GACAACCGGAATACCCTGCTGCTGTAGATGCTACTGCTGTGGGTGCCTGGCAGCACCATACTCTCTCGT" `shouldBe` 28.517294972946885
        score confWithoutForbidden "GACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCATACTCCCCCGTAACCGGAAC" `shouldBe` 59.71421227843634
        score confWithoutForbidden "GCCAGCGGCGACCCAAGACCCACACCTGTCCTAGCGGGT" `shouldBe` -12.47745413521671

        let forbidden = ["AACCGGAAC", "(A|C)..GGG(C|T)"]
        let confWithForbidden = CodonScoreConfig Human 3 1 1 0.5 1.4 40 0.001 2.6 23 1 43 forbidden
        score confWithForbidden "GACAACCGGAATACCCTGCTGCTGTAGATGCTACTGCTGTGGGTGCCTGGCAGCACCATACTCTCTCGT" `shouldBe` 28.517294972946885
        score confWithForbidden "GACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCATACTCCCCCGTAACCGGAAC" `shouldBe` 9.714212278436342
        score confWithForbidden "GCCAGCGGCGACCCAAGACCCACACCTGTCCTAGCGGGT" `shouldBe` -62.47745413521671

scoreFunDifferentForbiddenSeqWeight :: Spec
scoreFunDifferentForbiddenSeqWeight =
    describe "scoreFunDifferentForbiddenSeqWeight" $
    it "should correct count score with forbidden sequence weight = 4" $ do
        let conf' = CodonScoreConfig Human 3 1 1 0.5 1.4 40 0.001 2.6 100 4 43 defaultForbiddenRegexp
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCATACTCCCCCGGC" `shouldBe` -153.45210305442845
        score conf' "GCCAGCGGCGACCGGCGGGTCAAGACCCACACCTGTCCT" `shouldBe` -167.3461657496718
        score conf' "ACAGCCAGCGAATAAACCCCGAGGCCGCCGGCGGCCCTAGCGTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGACTACAGCCAAG" `shouldBe` -146.57214671578927
        score conf' "ACAGCCAGCGAATAAACCCCGACGATCGGGCCGCCGGCGAATGACAGCCAAG" `shouldBe` -137.875795016417

scoreFunDifferentGCDesired :: Spec
scoreFunDifferentGCDesired =
    describe "scoreFunDifferentGCDesired" $
    it "should correct count score with gc desired = 60" $ do
        let conf' = CodonScoreConfig Human 3 1 1 0.5 1.4 40 0.001 2.6 100 1 60 defaultForbiddenRegexp
        score conf' "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTGTGGGTGCCTGGCAGCACCGGC" `shouldBe` 80.83499151084192
        score conf' "GCCAGCGGCGACAAGACCCACACCTGTCCT" `shouldBe` 80.74362252733329
        score conf' "CCCTGCCCCGCCCCCGAGGCCGCCGGCGGCCCTAGCGTGTTCCTGTTCCCTCCTAAGCCTAAGGACACCCTGATGATCAGCAGAACCCCCGAGGTGACCTGCGTGGTGGTGGACGTGAGCCACGAGGACCCTGAGGTGAAGTTCAATTGGTACGTGGACGGCGTGGAGGTGCACAACGCCAAGACCAAGCCTAGAGAGGAGCAGTACAACAGCACCTACAGAGTGGTGAGCGTGCTGACCGTGCTGCACCAAGACTGGCTGAACGGCAAGGAGTACAAGTGCAAGGTGAGCAACAAGGCCCTGCCCGCCCCTATCGAGAAGACCATCAGCAAGGCCAAG"
         `shouldBe` 96.19662511761598
