module OligoDesigner.SpecGCContentOptimizer (
    gcContentOptimizerSpec
) where

import Test.Hspec (Spec, describe, it, shouldBe, shouldSatisfy)
import Bio.Tools.Sequence.OligoDesigner.Types (OligSplitting(..), OligSet(..), Olig(..), OligsDesignerInnerConfig(..))
import Bio.Tools.Sequence.OligoDesigner.Optimizer.GCContentOptimizer (gcContentOptimize)
import System.Random (mkStdGen)
import Bio.Tools.Sequence.CodonOptimization.Types (Organism(..))
import Bio.Tools.Sequence.OligoDesigner.Scorer (oligsGCContentDifference)
import Control.Monad.State.Lazy (evalStateT)
import Control.Monad.Except (runExcept)

gcContentOptimizerSpec :: Spec
gcContentOptimizerSpec =
    describe "gcContentOptimizerSpec" $ do
        optimizeGCContentSpec
        optimizeGCContentForDifferentTargetSpec
        optimizeGCContentRealDataSpec

optimizeGCContentSpec :: Spec
optimizeGCContentSpec =
    describe "optimizeGCContentSpec" $
    it "" $ do
        let oligs = OligSet
                        [Olig "TTGATCTTCC" 0 10, Olig "TTATAAGAAA" 10 20] -- 40% & 10%
                        [Olig "TATAAGGAAG" 5 15, Olig "TTGTTTTCT" 15 24]  -- 30% & 22%
                        (OligSplitting [(0, 10), (10, 20)] [(5, 15), (15, 24)])

        let gen = mkStdGen 6637
        let conf = OligsDesignerInnerConfig CHO 43 [] 0 0
        let (Right res) = runExcept $ evalStateT (gcContentOptimize conf oligs) gen
        oligsGCContentDifference res `shouldSatisfy` (<= oligsGCContentDifference oligs)
        res `shouldBe` OligSet
                        [Olig "TTGATCTTCC" 0 10, Olig "TTATACGGAA" 10 20] -- 40% & 30%
                        [Olig "TATAAGGAAG" 5 15, Olig "TTGTTTCCG" 15 24]  -- 30% & 40%
                        (OligSplitting [(0, 10), (10, 20)] [(5, 15), (15, 24)])

optimizeGCContentForDifferentTargetSpec :: Spec
optimizeGCContentForDifferentTargetSpec =
    describe "optimizeGCContentForDifferentTargetSpec" $
    it "" $ do
        let oligs = OligSet
                        [Olig "TTGATCTTCC" 0 10, Olig "TTATAAGAAA" 10 20] -- 40% & 10%
                        [Olig "TATAAGGAAG" 5 15, Olig "TTGTTTTCT" 15 24]  -- 30% & 22%
                        (OligSplitting [(0, 10), (10, 20)] [(5, 15), (15, 24)])

        let gen = mkStdGen 6637
        let conf = OligsDesignerInnerConfig CHO 20 [] 0 0
        let (Right res) = runExcept $ evalStateT (gcContentOptimize conf oligs) gen
        oligsGCContentDifference res `shouldSatisfy` (<= oligsGCContentDifference oligs)
        res `shouldBe` OligSet
                        [Olig "TTGATCTTCT" 0 10, Olig "TGATAAGAAA" 10 20] -- 30% & 20%
                        [Olig "TATCAAGAAG" 5 15, Olig "TTGTTTTCT" 15 24]  -- 30% & 22%
                        (OligSplitting [(0, 10), (10, 20)] [(5, 15), (15, 24)])

optimizeGCContentRealDataSpec :: Spec
optimizeGCContentRealDataSpec =
    describe "optimizeGCContentRealDataSpec" $
    it "" $ do
        let coords = OligSplitting [(0, 60), (60, 120), (120, 180), (180, 240)] [(30, 90), (90, 150), (150, 210), (210, 270)]
        let oligs =
                OligSet
                    [ Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60      -- 66%
                    , Olig "GGCACCGCCGCCCTGGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120    -- 68%
                    , Olig "TGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGC" 120 180   -- 70%
                    , Olig "GGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACC" 180 240   -- 68%
                    ]
                    [ Olig "CTTCACCAGGCAGCCCAGGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90     -- 73%
                    , Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150    -- 66%
                    , Olig "GGTCACCACGCTGCTCAGGCTGTACAGGCCGCTGCTCTGCAGCACGGCAGGGAAGGTGTG" 150 210   -- 66%
                    , Olig "GCTAGGCTTGTGGTTCACGTTGCAGATGTAGGTCTGGGTGCCCAGGCTGCTGCTAGGCAC" 210 270   -- 60%
                    ]
                    coords
        let gen = mkStdGen 6637
        let conf = OligsDesignerInnerConfig CHO 43 [] 0 0
        let (Right res) = runExcept $ evalStateT (gcContentOptimize conf oligs) gen
        oligsGCContentDifference res `shouldSatisfy` (<= oligsGCContentDifference oligs)
        res `shouldBe` OligSet
                    [ Olig "GCTAGCACCAAGGGCCCCAGCGTGTTTCCTCTGGCCCCTAGCAGCAAGAGCACCAGCGGC" 0 60      -- 66%
                    , Olig "GGCACCGCCGCCTTAGGCTGCCTGGTGAAGGACTACTTCCCTGAGCCTGTGACCGTGAGC" 60 120    -- 65%
                    , Olig "TGGAACAGCGGCGCCCTGACCAGCGGCGTGCACACCTTCCCTGCCGTGCTGCAGAGCAGC" 120 180   -- 70%
                    , Olig "GGCCTGTACAGCCTGAGCAGCGTGGTGACCGTGCCTAGCAGCAGCCTGGGCACCCAGACC" 180 240   -- 68%
                    ]
                    [ Olig "CTTCACCAGGCAGCCTAAGGCGGCGGTGCCGCCGCTGGTGCTCTTGCTGCTAGGGGCCAG" 30 90     -- 70%
                    , Olig "CACGCCGCTGGTCAGGGCGCCGCTGTTCCAGCTCACGGTCACAGGCTCAGGGAAGTAGTC" 90 150    -- 66%
                    , Olig "GGTCACCACGCTGCTCAGGCTGTACAGGCCGCTGCTCTGCAGCACGGCAGGGAAGGTGTG" 150 210   -- 66%
                    , Olig "GCTAGGCTTGTGGTTCACGTTGCAGATGTAGGTCTGGGTGCCCAGGCTGCTGCTAGGCAC" 210 270   -- 60%
                    ]
                    coords

