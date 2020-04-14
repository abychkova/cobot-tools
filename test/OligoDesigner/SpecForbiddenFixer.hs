module OligoDesigner.SpecForbiddenFixer (
    fixerSpec
) where

import Control.Monad.Except     (runExcept)
import Control.Monad.State.Lazy (evalStateT)
import System.Random            (mkStdGen)
import Test.Hspec               (Spec, describe, it, shouldBe)
import Text.Regex.TDFA          (Regex, makeRegex, match)

import Bio.Tools.Sequence.CodonOptimization.Types        (Organism (..))
import Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer   (fixForbidden)
import Bio.Tools.Sequence.OligoDesigner.Types            (OligsDesignerInnerConfig (..))
import Bio.Tools.Sequence.OligoDesigner.Utils.Prettifier (prettyDNA)

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
        let conf = OligsDesignerInnerConfig CHO 60 [regexp] 0 5
        let dna = "GCCAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTTCTTCTAAGTCTACCTCTGGCGGCACCGCCGCCCTGGGCTGTCTGGTGAAG"
        
        let (Right res) = runExcept $ evalStateT (fixForbidden conf dna) gen
        let isMatch = regexp `match` prettyDNA res :: Bool
        isMatch `shouldBe` False
        
forbiddenConstantFixerSpec :: Spec
forbiddenConstantFixerSpec =
    describe "forbiddenConstantFixerSpec" $
    it "" $ do
        let gen = mkStdGen 429
        let regexpStr = "TGG" :: String
        let regexp = makeRegex regexpStr :: Regex
        let conf = OligsDesignerInnerConfig CHO 60 [regexp] 0 5
        let dna = "GCCTGGACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTTCTTCTAAGTCTACCTCTGGCGGCACCGCCGCCCTGGGCTGTCTGGTGAAG"
        
        let (Left msg) = runExcept $ evalStateT (fixForbidden conf dna) gen
        msg `shouldBe` "cannot fix this shit"
        
fixTowForbiddenSpec :: Spec
fixTowForbiddenSpec =
    describe "fixTowForbiddenSpec" $
    it "" $ do
        let gen = mkStdGen 3359
        let regexpStr = ["CAGG", "AATAAA", "GCCGCCATGG"] :: [String]
        let regexps = map makeRegex regexpStr :: [Regex]
        let conf = OligsDesignerInnerConfig CHO 60 regexps 0 5
        let dna = "CAGGCCGCCATGGGCAATAAACAGGTG" --CAGG  AATAAA  GCCGCCATGG

        let (Right res) = runExcept $ evalStateT (fixForbidden conf dna) gen
        let matches = [regexp `match` prettyDNA res :: Bool | regexp <- regexps]
        filter (== True) matches `shouldBe` []
        
fixTowForbiddenRealDataSpec :: Spec
fixTowForbiddenRealDataSpec =
    describe "fixTowForbiddenRealDataSpec" $
    it "" $ do
        let gen = mkStdGen 3358
        let regexpStr = ["CAGG", "AATAAA", "GCCGCCATGG"] :: [String]
        let regexps = map makeRegex regexpStr :: [Regex]
        let conf = OligsDesignerInnerConfig CHO 60 regexps 0 5
        let dna = "CAGGTGCAGCTGCAGGAGAGCGGCGGCGGCCTGGTGCAGCCTGGCGGCTCTCTGAGACTGTCTTGTGCCGCCTCTGGCATCCAGTTCAAGTTCCACAACATGGGCTGGGGCAGACAAGCCCCTGGCAAGCAGAGAGAGCACGTGGCCGCCATGGATCACGGCAGAAGAACCATCTACGCCGACCACGTGAAGGGCAGATTCACCATCTCTAGAGACAACACCAGGAACACCGTGTACCTGCAGATGAACTCTCTGAAGGCCGAGGACACCGCCATCTACTACTGCAAGGCCTCTGCCGGCAGAAGAGTGTACTGGGGCCAAGGCACCATGGTGACCGTGTCTTCT" --CAGG  AATAAA  GCCGCCATGG

        let (Right res) = runExcept $ evalStateT (fixForbidden conf dna) gen
        let matches = [regexp `match` prettyDNA res :: Bool | regexp <- regexps]
        filter (== True) matches `shouldBe` []
