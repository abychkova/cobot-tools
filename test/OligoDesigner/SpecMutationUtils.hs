module OligoDesigner.SpecMutationUtils (
    mutationUtilsSpec
) where

import Bio.Tools.Sequence.OligoDesigner.Utils.MutationUtils (weightedRandom, randomCodon, oneMutation, mutateSlice, mutate)
import           Test.Hspec                             (Spec, describe,
                                                         errorCall, it,
                                                         shouldBe,
                                                         shouldSatisfy,
                                                         shouldThrow, shouldNotBe)
import System.Random (StdGen, getStdGen, mkStdGen)
import Control.Monad.State.Lazy (State, get)
import Bio.Tools.Sequence.CodonOptimization.Types (Organism(..))
import Control.Monad.State (evalState, runState, put)
import Bio.NucleicAcid.Nucleotide.Type (DNA(..))
import Control.Exception (evaluate)
import Bio.Protein.AminoAcid (AA(..))
import Data.List (nub)

mutationUtilsSpec :: Spec
mutationUtilsSpec =
    describe "mutationUtilsSpec" $ do
        randomCodonSpec
        randomCodonForAAWithOneCodonSpec
        randomCodonWeightedSpec
        randomCodonForDifferentOrganismSpec
        
        oneMutationSpec
        oneMutationForAAWithOneCodonSpec

        mutateSliceSpec
        mutateSliceWhenThereIsNoVariantsSpec
        mutateSliceRealRandomSpec

        mutateSpec
        mutateFromStartSpec
        mutateInvalidIntervalSpec
        mutateOneAASpec

        weightedRandomSimpleSpec
        weightedRandomMultipleCallSpec
        weightedRandomSpec
        weightedRandomAllZeroSpec
        weightedRandomForUnsortedSpec
        weightedRandomForEmptySpec
        
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
        elem "TCTTTGCCGAACGAGGGCATG" res `shouldBe` True
        elem "AGCTTGCCGAACGAGGGCATG" res `shouldBe` True
        elem "TCTCTCCCGAACGAGGGCATG" res `shouldBe` True
        elem "TCTTTGCCAAACGAGGGCATG" res `shouldBe` True
        elem "TCTTTGCCGAATGAGGGCATG" res `shouldBe` True
        elem "TCTTTGCCGAACGAAGGCATG" res `shouldBe` True
        elem "TCTTTGCCGAACGAGGGAATG" res `shouldBe` True

--TODO: test me with forbidden regexp
mutateSpec :: Spec
mutateSpec =
    describe "mutateSpec" $
    it "" $ do
        let dna = "ATGGAGACC" ++ "AATATGCAT" ++ "GACACCCTGCTGCTGTGGGTGCTGCTGCTG"
        let gen = mkStdGen 4
        let res = evalState (mutate CHO dna (4, 6)) gen
        res `shouldBe` ["ATGGAGACCAATATGCATGACACCCTGCTGCTGTGGGTGCTGCTGCTG",
                        "ATGGAGACCAACATGCATGACACCCTGCTGCTGTGGGTGCTGCTGCTG",
                        "ATGGAGACCAATATGCACGACACCCTGCTGCTGTGGGTGCTGCTGCTG"]

mutateFromStartSpec :: Spec
mutateFromStartSpec =
    describe "mutateFromStartSpec" $
    it "" $ do
        let dna = "AATATGCAT" ++ "ATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTG"
        let gen = mkStdGen 4
        let res = evalState (mutate CHO dna (1, 3)) gen
        res `shouldBe` ["AATATGCATATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTG",
                        "AACATGCATATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTG",
                        "AATATGCACATGGAGACCGACACCCTGCTGCTGTGGGTGCTGCTGCTG"]

mutateOneAASpec :: Spec
mutateOneAASpec =
    describe "mutateOneAASpec" $
    it "" $ do
        let dna = "ATGGAGACC" ++ "AAT" ++ "ATGCATGACACCCTGCTGCTGTGGGTGCTGCTGCTG"
        let gen = mkStdGen 4
        let res = evalState (mutate CHO dna (4, 4)) gen
        res `shouldBe` ["ATGGAGACCAATATGCATGACACCCTGCTGCTGTGGGTGCTGCTGCTG",
                        "ATGGAGACCAACATGCATGACACCCTGCTGCTGTGGGTGCTGCTGCTG"]

mutateInvalidIntervalSpec :: Spec
mutateInvalidIntervalSpec =
    describe "mutateInvalidIntervalSpec" $
    it "" $ do
        let gen = mkStdGen 4
        evaluate (evalState (mutate CHO "AATATGCATATG" (3, 1)) gen) `shouldThrow` errorCall "invalid interval for mutation: (3,1)"
        evaluate (evalState (mutate CHO "AATATGCATATG" (-1, 1)) gen) `shouldThrow` errorCall "invalid interval for mutation: (-1,1)"
        evaluate (evalState (mutate CHO "AATATGCATATG" (3, -10)) gen) `shouldThrow` errorCall "invalid interval for mutation: (3,-10)"
        evaluate (evalState (mutate CHO "AATATGCATATG" (3, 5)) gen) `shouldThrow` errorCall "invalid interval for mutation: (3,5)"
        evaluate (evalState (mutate CHO "AATATGCATATG" (5, 40)) gen) `shouldThrow` errorCall "invalid interval for mutation: (5,40)"

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
