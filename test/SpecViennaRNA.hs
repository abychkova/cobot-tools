module SpecViennaRNA where

import           Bio.NucleicAcid.Nucleotide          (DNA)
import           Bio.Tools.Sequence.ViennaRNA.Cofold (cofold)
import           Bio.Tools.Sequence.ViennaRNA.Fold   (fold)
import           Test.Hspec

foldSpec :: Spec
foldSpec = it "should work" $ do
    let (energy1, structure1) = fold 37 ("CTGGATCGCAATGACGCTCTTAGGTCTCGT" :: [DNA])
    energy1 `shouldBe` -2.9
    structure1 `shouldBe` "..(((((((......)).....)))))..."

    let (energy2, structure2) = fold 37 ("CTGGATCGCAATGGGTCTCGT" :: [DNA])
    energy2 `shouldBe` -4.4
    structure2 `shouldBe` "..(((((......)))))..."

cofoldSpec :: Spec
cofoldSpec = it "should work" $ do
    let (energy1, structure1) = cofold 37 ("CAAGTACAGTTACAAGAAAGTGGAGGAGGATTAGTACAACCGGGAGGAAGTCTCAG" :: [DNA], "TCTGAATCCACTCGCAGCACAGGAGAGTCTGAGACTTCCTCCCGGTTGTACTAATC" :: [DNA])
    energy1 `shouldBe` -56.40
    structure1 `shouldBe` "......((.((......)).))......((((((((((((((((((((((((((((.........((((...........))))))))))))))))))))))))))))))))"

    let (energy2, structure2) = cofold 37 ("TCTGAATCCACTCGCAGCACAGGAGAGTCTGAGACTTCCTCCCGGTTGTACTAATC" :: [DNA], "ACTCTCCTGTGCTGCGAGTGGATTCAGATTCAGTAACTACGCGATGAGTTGGGTCC" :: [DNA])
    energy2 `shouldBe` -62.30
    structure2 `shouldBe` "((((((((((((((((((((((((((((...((((..((....))..)).))....))))))))))))))))))))))))))))..........((.(((....))).)).."

foldTest :: Spec
foldTest = describe "Predict RNA secondary structure and calculate energy for one sequence" $ do
    foldSpec

cofoldTest :: Spec
cofoldTest = describe "Predict RNA secondary structure and calculate energy for two sequences" $ do
    cofoldSpec
