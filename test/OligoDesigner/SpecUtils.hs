{-# LANGUAGE OverloadedStrings #-}

module OligoDesigner.SpecUtils where

import           Bio.Tools.Sequence.OligoDesigner.Types (Olig (..),
                                                         OligSet (..),
                                                         OligSplitting (..))
import           Bio.Tools.Sequence.OligoDesigner.Utils (assemble,
                                                         weightedRandom)
import           Control.Exception                      (evaluate)
import           Control.Monad.State                    (State, evalState, get,
                                                         put, runState)
import           Data.List                              (nub)
import           System.Random                          (StdGen, getStdGen)
import           Test.Hspec                             (Spec, describe,
                                                         errorCall, it,
                                                         shouldBe,
                                                         shouldSatisfy,
                                                         shouldThrow)

utilsSpec :: Spec
utilsSpec =
    describe "utilsSpec" $ do
        assembleSpec
        assembleWithGapSpec
        weightedRandomSimpleSpec
        weightedRandomMultipleCallSpec
        weightedRandomSpec
        weightedRandomAllZeroSpec
        weightedRandomForUnsortedSpec
        weightedRandomForEmptySpec

assembleSpec :: Spec
assembleSpec =
    describe "assembleSpec" $
    it "assembleSpec" $ do
        let fwd = [Olig "AA" 0 2, Olig "TT" 2 4, Olig "GG" 4 6, Olig "CC" 6 8, Olig "AA" 8 10]
        let rsd = [Olig "TA" 1 3, Olig "AC" 3 5, Olig "CG" 5 7, Olig "GT" 7 9, Olig "TA" 9 11]
        let oligs = OligSet fwd rsd (OligSplitting [] [])
        assemble oligs `shouldBe` "AATTGGCCAAT"

assembleWithGapSpec :: Spec
assembleWithGapSpec =
    describe "assembleWithGapSpec" $
    it "assembleWithGapSpec" $ do
        let fwd = [Olig "AATT" 0 4, Olig "CCAA" 6 10]
        let rsd = [Olig "ACCG" 3 7, Olig "TTA" 8 11]
        let oligs = OligSet fwd rsd (OligSplitting [] [])
        assemble oligs `shouldBe` "AATTGGCCAAT"

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
