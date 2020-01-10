{-# LANGUAGE OverloadedStrings #-}

module OligoDesigner.SpecUtils where

import Bio.NucleicAcid.Nucleotide (DNA)
import Bio.Tools.Sequence.OligoDesigner.Types (Olig(..), OligSet(..))
import Bio.Tools.Sequence.OligoDesigner.Utils (assemble)
import Bio.Tools.Sequence.OligoDesigner.Algo (generateOligs)
import Test.Hspec (Spec, describe, it, shouldBe)
import Debug.Trace (trace)

utilsSpec :: Spec
utilsSpec = describe "utilsSpec" $ do
    assembleSpec
    assembleWithGapSpec

assembleSpec :: Spec
assembleSpec =
    describe "assembleSpec" $
    it "assembleSpec" $
     do
        let fwd = [Olig "AA" 0 2, Olig "TT" 2 4, Olig "GG" 4 6, Olig "CC" 6 8, Olig "AA" 8 10]
        let rsd = [Olig "TA" 1 3, Olig "AC" 3 5, Olig "CG" 5 7, Olig "GT" 7 9, Olig "TA" 9 11]
        let oligs = OligSet fwd rsd
        assemble oligs `shouldBe` "AATTGGCCAAT"

assembleWithGapSpec :: Spec
assembleWithGapSpec =
    describe "assembleWithGapSpec" $
    it "assembleWithGapSpec" $
     do
        let fwd = [Olig "AATT" 0 4, Olig "CCAA" 6 10]
        let rsd = [Olig "ACCG" 3 7, Olig "TTA" 8 11]
        let oligs = OligSet fwd rsd
        assemble oligs `shouldBe` "AATTGGCCAAT"