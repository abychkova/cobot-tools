module OligoDesigner.SpecIterationOptimizer (
    optimizerSpec
) where

import Test.Hspec (Spec, describe, it, shouldBe, shouldNotBe, shouldSatisfy)
import Control.Monad.State (evalState)
import System.Random (mkStdGen)
import Bio.Tools.Sequence.OligoDesigner.Optimizer.IterativeOptimizer (optimize)
import Bio.Tools.Sequence.OligoDesigner.Scorer (commonScore)
import Bio.Tools.Sequence.OligoDesigner.Types (OligsDesignerConfig(..), OligsSplittingConfig(..), OligSet(..), Olig(..),
    OligSplitting(..), OligsDesignerInnerConfig(..))
import Data.Default (def)
import Bio.Tools.Sequence.CodonOptimization.Types (Organism(..))

--TODO: test me with forbidden regexp
optimizerSpec :: Spec
optimizerSpec =
    describe "optimizerSpec" $ do
        optimizeWithZeroIterationsSpec
        optimizeWithOneIterationsSpec
        optimizeTillStableScoreSpec

optimizeWithZeroIterationsSpec :: Spec
optimizeWithZeroIterationsSpec =
    describe "optimizeWithZeroIterationsSpec" $
    it "" $ do
        let conf = OligsDesignerInnerConfig CHO 43 [] 0 1
        let gen = mkStdGen 63
        let res = evalState (optimize conf oligs) gen
        res `shouldBe` oligs

optimizeWithOneIterationsSpec :: Spec
optimizeWithOneIterationsSpec =
    describe "optimizeWithOneIterationsSpec" $
    it "" $ do
        let conf = OligsDesignerInnerConfig CHO 43 [] 1 1
        let gen = mkStdGen 63
        let res = evalState (optimize conf oligs) gen
        res `shouldNotBe` oligs
        commonScore 43 res `shouldSatisfy` (> commonScore 43 oligs)

optimizeTillStableScoreSpec :: Spec
optimizeTillStableScoreSpec =
    describe "optimizeTillStableScoreSpec" $
    it "" $ do
        let conf = OligsDesignerInnerConfig CHO 43 [] 8 1
        let gen = mkStdGen 63
        let res = evalState (optimize conf oligs) gen
        res `shouldNotBe` oligs
        commonScore 43 res `shouldSatisfy` (> commonScore 43 oligs)

oligs :: OligSet
oligs = OligSet [ Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60
                , Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120
                , Olig "TGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGC" 120 180
                , Olig "GGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACC" 180 240
                ]
                [ Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90
                , Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150
                , Olig "GGTCACCACGCTGCTCAGGCTGTACAGGCCGCTGCTCTGCAGCACGGCAGGGAAGGTGTG" 150 210
                , Olig "GCTAGGCTTGTGGTTCACGTTGCAGATGTAGGTCTGGGTGCCCAGGCTGCTGCTAGGCAC" 210 270
                ]
                (OligSplitting [(0, 60), (60, 120), (120, 180), (180, 240)] [(30, 90), (90, 150), (150, 210), (210, 270)])