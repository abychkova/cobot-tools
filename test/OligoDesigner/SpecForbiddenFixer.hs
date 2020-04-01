{-# LANGUAGE OverloadedStrings #-}

module OligoDesigner.SpecForbiddenFixer where

import Test.Hspec (Spec, describe, it, shouldBe)
import Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer (fixForbidden)
import Text.Regex.TDFA (match, makeRegex, Regex, (=~))
import Bio.Tools.Sequence.OligoDesigner.Prettifier (prettyDNA)
import System.Random (mkStdGen)
import Bio.Tools.Sequence.CodonOptimization (CodonOptimizationConfig(..))
import Bio.Tools.Sequence.CodonOptimization.Types (Organism(..))
import Bio.Tools.Sequence.OligoDesigner.Types (OligsSplittingConfig(..), OligsDesignerConfig(..))
import Data.Default (def)
import Control.Monad.Except (runExcept)
import Debug.Trace (trace)

fixerSpec :: Spec
fixerSpec =
    describe "fixerSpec" $ do
        forbiddenFixerSpec
        forbiddenConstantFixerSpec
        fixTowForbiddenSpec
        fixTowForbiddenRealDataSpec

forbiddenFixerSpec :: Spec
forbiddenFixerSpec =
    describe "forbiddenFixerSpec" $
    it "" $ do
        let gen = mkStdGen 429
        let regexpStr = "GGGCCCC" :: String
        let regexp = makeRegex regexpStr :: Regex
        let conf = OligsDesignerConfig def def 0 0 0 0 5
        let dna = "GCCAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTTCTTCTAAGTCTACCTCTGGCGGCACCGCCGCCCTGGGCTGTCTGGTGAAG"
        
        let (Right res) = runExcept $ fixForbidden gen conf [regexp] dna
        let isMatch = regexp `match` prettyDNA res :: Bool
        trace (prettyDNA res) $ isMatch `shouldBe` False
        
forbiddenConstantFixerSpec :: Spec
forbiddenConstantFixerSpec =
    describe "forbiddenConstantFixerSpec" $
    it "" $ do
        let gen = mkStdGen 429
        let regexpStr = "TGG" :: String
        let regexp = makeRegex regexpStr :: Regex
        let conf = OligsDesignerConfig def def 0 0 0 0 5
        let dna = "GCCTGGACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTTCTTCTAAGTCTACCTCTGGCGGCACCGCCGCCCTGGGCTGTCTGGTGAAG"
        
        let (Left msg) = runExcept $ fixForbidden gen conf [regexp] dna
        msg `shouldBe` "cannot fix this shit"
        
fixTowForbiddenSpec :: Spec
fixTowForbiddenSpec =
    describe "fixTowForbiddenSpec" $
    it "" $ do
        let gen = mkStdGen 3359
        let regexpStr = ["CAGG", "AATAAA", "GCCGCCATGG"] :: [String]
        let regexps = map makeRegex regexpStr :: [Regex]
        let conf = OligsDesignerConfig def def 0 0 0 0 5
        let dna = "CAGGCCGCCATGGGCAATAAACAGGTG" --CAGG  AATAAA  GCCGCCATGG

        let (Right res) = runExcept $ fixForbidden gen conf regexps dna
        let matches = [regexp `match` prettyDNA res :: Bool | regexp <- regexps]
        trace (prettyDNA res) $ filter (== True) matches `shouldBe` []
        
fixTowForbiddenRealDataSpec :: Spec
fixTowForbiddenRealDataSpec =
    describe "fixTowForbiddenRealDataSpec" $
    it "" $ do
        let gen = mkStdGen 3358
        let regexpStr = ["CAGG", "AATAAA", "GCCGCCATGG"] :: [String]
        let regexps = map makeRegex regexpStr :: [Regex]
        let conf = OligsDesignerConfig def def 0 0 0 0 5
        let dna = "CAGGTGCAGCTGCAGGAGAGCGGCGGCGGCCTGGTGCAGCCTGGCGGCTCTCTGAGACTGTCTTGTGCCGCCTCTGGCATCCAGTTCAAGTTCCACAACATGGGCTGGGGCAGACAAGCCCCTGGCAAGCAGAGAGAGCACGTGGCCGCCATGGATCACGGCAGAAGAACCATCTACGCCGACCACGTGAAGGGCAGATTCACCATCTCTAGAGACAACACCAGGAACACCGTGTACCTGCAGATGAACTCTCTGAAGGCCGAGGACACCGCCATCTACTACTGCAAGGCCTCTGCCGGCAGAAGAGTGTACTGGGGCCAAGGCACCATGGTGACCGTGTCTTCT" --CAGG  AATAAA  GCCGCCATGG

        let (Right res) = runExcept $ fixForbidden gen conf regexps dna
        let matches = [regexp `match` prettyDNA res :: Bool | regexp <- regexps]
        trace (prettyDNA res) $ filter (== True) matches `shouldBe` []