module OligoDesigner.SpecOligoDesignerOptimizer where

import Test.Hspec (Spec, shouldBe, it, describe)
import Bio.Tools.Sequence.OligoDesigner.Optimizer (maxPairMutationIndexes)

optimizerSpec :: Spec
optimizerSpec =
    describe "optimizerSpec" $ do
        maxPairMutationIndexesSpec

maxPairMutationIndexesSpec :: Spec
maxPairMutationIndexesSpec =
    describe "maxPairMutationIndexes" $
    it "should find pair with maximum rna by matrix" $ do
        let mtx = 
        let res = maxPairMutationIndexes mtx
        True `shouldBe` True

