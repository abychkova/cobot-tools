module SpecPrimers where

import           Bio.NucleicAcid.Nucleotide              (DNA, symbol)
import           Bio.Tools.Sequence.Primers.Optimization
import           Bio.Tools.Sequence.Primers.Types
import           Test.Hspec

testPrimers :: SpecWith ()
testPrimers = describe "Primer optimization test" testPrimersOptimization

testPrimersOptimization :: SpecWith ()
testPrimersOptimization = do
    it "Sequence: GATCCAACTTCAAAGAGTCCTGGCCAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGC; pos: 0" $
      let s = "GATCCAACTTCAAAGAGTCCTGGCCAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGC" :: [DNA]
      in fmap symbol . _seq' <$> designPrimer False s 0 `shouldBe` Right "GATCCAACTTCAAAGAGTCCTGGC"
    it "Sequence: GATCCAACTTCAAAGAGTCCTGGCCAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGC; pos: 61" $
      let s = "GATCCAACTTCAAAGAGTCCTGGCCAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGC" :: [DNA]
      in fmap symbol . _seq' <$> designPrimer False s 61 `shouldBe` Right "CACAACTAGAATGCAGTGAAAAAAATG"
    it "Sequence: GATCCAACTTCAAAGAGTCCTGGCCAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGC; pos: 90" $
      let s = "GATCCAACTTCAAAGAGTCCTGGCCAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGC" :: [DNA]
      in fmap symbol . _seq' <$> designPrimer False s 90 `shouldBe` Right "TTATTTGTGAAATTTGTGATGC"
    it "Sequence: GATCCAACTTCAAAGAGTCCTGGCCAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGCTATTGCTTTATTTGTAACCATTATAAGCTGCAATAAACAAGTTTCAGGCACCGGGCTTGCGGGTCATGCAC; pos: 90" $
      let s = "GATCCAACTTCAAAGAGTCCTGGCCAGACATGATAAGATACATTGATGAGTTTGGACAAACCACAACTAGAATGCAGTGAAAAAAATGCTTTATTTGTGAAATTTGTGATGCTATTGCTTTATTTGTAACCATTATAAGCTGCAATAAACAAGTTTCAGGCACCGGGCTTGCGGGTCATGCAC" :: [DNA]
      in fmap symbol . _seq' <$> designPrimer False s 90 `shouldBe` Right "TTATTTGTGAAATTTGTGATGCTATTGC"
    it "Shouldn't find any primers. Sequence: GCAATAGCATCACAAATTTCACAAATAAGCAATAGCATCACAAATTTCACAAATAAGCAATAGCATCACAAATTTCACAAATAATTATTTGTGAAATTTGTGATGCTATTGCTTATTTGTGAAATTTGTGATGCTATTGCTTATTTGTGAAATTTGTGATGCTATTGC; pos: 0" $
      let s = "GCAATAGCATCACAAATTTCACAAATAAGCAATAGCATCACAAATTTCACAAATAAGCAATAGCATCACAAATTTCACAAATAATTATTTGTGAAATTTGTGATGCTATTGCTTATTTGTGAAATTTGTGATGCTATTGCTTATTTGTGAAATTTGTGATGCTATTGC" :: [DNA]
      in fmap symbol . _seq' <$> designPrimer False s 0 `shouldBe` Left "Bio.Tools.Sequence.Primers.Optimization: all primers designed from given position don't bind to target."
    it "Shouldn't find any primers. Sequence: ACAAATAATTATTTGTGAAATTTGACAAATAATTATTTGTGAAATTTGACAAATAATTATTTGTGAAATTTGCAAATTTCACAAATAATTATTTGTCAAATTTCACAAATAATTATTTGTCAAATTTCACAAATAATTATTTGT; pos: 61" $
      let s = "GCAATAGCATCACAAATTTCACAAATAAGCAATAGCATCACAAATTTCACAAATAAGCAATAGCATCACAAATTTCACAAATAATTATTTGTGAAATTTGTGATGCTATTGCTTATTTGTGAAATTTGTGATGCTATTGCTTATTTGTGAAATTTGTGATGCTATTGC" :: [DNA]
      in fmap symbol . _seq' <$> designPrimer False s 61 `shouldBe` Left "Bio.Tools.Sequence.Primers.Optimization: all primers designed from given position don't bind to target."
    let wholePlasmid = "CCTGCAGGCAGCTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCTGCGGCCCGACATTCCGGAGGTACCTCTAGATTGGCAAACAGCTATTATGGGTATTATGGGTGATCTCCAGATGGCTAAACTTTTAAATCATGAATGAAGTAGATATTACCAAATTGCTTTTTCAGCATCCATTTAGATAATCATGTTTTTTGCCTTTAATCTGTTAATGTAGTGAATTACAGAAATACATTTCCTAAATCATTACATCCCCCAAATCGTTAATCTGCTAAAGTACATCTCTGGCTCAAACAAGACTGGTTGTGCATCTCAATTAGTCAGCAACCATAGTCCCGCCCCTAACTCCGCCCATCCCGCCCCTAACTCCGCCCAGTTCCGCCCATTCTCCGCCCCATCGCTGACTAATTTTTTTTATTTATGCAGAGGCCGAGGCCGCCTCGGCCTCTGAGCTATTCCAGAAGTAGTGAGGAGGCTTTTTTGGAGGCCTAGGCTTTTGCAAACGTAACTATAACGGTCCTAAGGTAGCGAAATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGACCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTGATAGGGATAACAGGGTAATGTCGACAAGCTTGCTAGCACGGGTGGCATCCCTGTGACCCCTCCCCAGTGCCTCTCCTGGCCCTGGAAGTTGCCACTCCAGTGCCCACCAGCCTTGTCCTAATAAAATTAAGTTGCATCATTTTGTCTGACTAGGTGTCCTTCTATAATATTATGGGGTGGAGGGGGGTGGTATGGAGCAAGGGGCAAGTTGGGAAGACAACCTGTAGGGCCTGCGGGGTCTATTGGGAACCAAGCTGGAGTGCAGTGGCACAATCTTGGCTCACTGCAATCTCCGCCTCCTGGGTTCAAGCGATTCTCCTGCCTCAGCCTCCCGAGTTGTTGGGATTCCAGGCATGCATGACCAGGCTCAGCTAATTTTTGTTTTTTTGGTAGAGACGGGGTTTCACCATATTGGCCAGGCTGGTCTCCAACTCCTAATCTCAGGTGATCTACCCACCTTGGCCTCCCAAATTGCTGGGATTACAGGCGTGAACCACTGCTCCCTTCCCTGTCCTTATGCATGGGCCCGTACGATCACTAGTGTACAGCGGCCGCAGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAGCTGCCTGCAGGGGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATACGTCAAAGCAACCATAGTACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAGCGCCTTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTCTCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTAGGGTTCCGATTTAGTGCTTTACGGCACCTCGACCCCAAAAAACTTGATTTGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATAGACGGTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACTCTATCTCGGGCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGTCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAACGCGAATTTTAACAAAATATTAACGTTTACAATTTTATGGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAGTTAAGCCAGCCCCGACACCCGCCAACACCCGCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGGCATCCGCTTACAGACAAGCTGTGACCGTCTCCGGGAGCTGCATGTGTCAGAGGTTTTCACCGTCATCACCGAAACGCGCGAGACGAAAGGGCCTCGTGATACGCCTATTTTTATAGGTTAATGTCATGATAATAATGGTTTCTTAGACGTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAGACAATAACCCTGATAAATGCTTCAATAATATTGAAAAAGGAAGAGTATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGT"
    it "Shouldn't find any primers; Sequence: whole plasmid; pos: 5" $
      fmap symbol . _seq' <$> designPrimer True wholePlasmid 5 `shouldBe` Left "Bio.Tools.Sequence.Primers.Optimization: all primers designed from given position don't bind to target."
    it "Sequence: whole plasmid; pos: 1781" $
      fmap symbol . _seq' <$> designPrimer True wholePlasmid 1781 `shouldBe` Right "TGATCTACCCACCTTGGCCTCCC"
    it "Sequence: whole plasmid; pos: 4628 (last index in linear sequence)" $
      fmap symbol . _seq' <$> designPrimer True wholePlasmid 4628 `shouldBe` Right "TCCTGCAGGCAGCTGCGCGC"
    it "Should fail, because index is out of range; Sequence: whole plasmid; pos: 5000" $
      fmap symbol . _seq' <$> designPrimer True wholePlasmid 5000 `shouldBe` Left "Bio.Tools.Sequence.Primers.Optimization: given position is out of range."
    let wholePlasmid2 = "GGACGTCCGTCGACGCGCGAGCGAGCGAGTGACTCCGGCGGGCCCGCAGCCCGCTGGAAACCAGCGGGCCGGAGTCACTCGCTCGCTCGCGCGTCTCTCCCTCACCGGTTGAGGTAGTGATCCCCAAGGATGCGCAGATCAATAATTATCATTAGTTAATGCCCCAGTAATCAAGTATCGGGTATATACCTCAAGGCGCAATGTATTGAATGCCATTTACCGGGCGGACCGACTGGCGGGTTGCTGGGGGCGGGTAACTGCAGTTATTACTGCATACAAGGGTATCATTGCGGTTATCCCTGAAAGGTAACTGCAGTTACCCACCTCATAAATGCCATTTGACGGGTGAACCGTCATGTAGTTCACATAGTATACGGTTCATGCGGGGGATAACTGCAGTTACTGCCATTTACCGGGCGGACCGTAATACGGGTCATGTACTGGAATACCCTGAAAGGATGAACCGTCATGTAGATGCATAATCAGTAGCGATAATGGTACCACTACGCCAAAACCGTCATGTAGTTACCCGCACCTATCGCCAAACTGAGTGCCCCTAAAGGTTCAGAGGTGGGGTAACTGCAGTTACCCTCAAACAAAACCGTGGTTTTAGTTGCCCTGAAAGGTTTTACAGCATTGTTGAGGCGGGGTAACTGCGTTTACCCGCCATCCGCACATGCCACCCTCCAGATATATTCGTCTCGAGCAAATCACTTGGCAGTCTAGCGGACCTCTGCGGTAGGTGCGACAAAACTGGAGGTATCTTCTGTGGCCCTGGCTAGGTCGGAGGTACGTCGCGCACTTGTACTAGTACCGGCTCTCGGGACCGGACTAGTGGTAGACGGACGACCCGATGGACGACTCGCGGCTCACGTGGCACAAGGACCTGGTGCTCTTGCGGTTGTTCTAGGACTTGGCCGGGTTCTCTATGTTGTCGCCGTTCGACCTCCTCAAGCACGTCCCGTTGGACCTCTCCCTCACGTACCTCCTCTTCACGTCGAAGCTCCTCCGGTCCCTTCACAAGCTCTTGTGGCTCGCCTGGTGGCTCAAGACCTTCGTCATGCACCTGCCGCTGGTCACGCTCTCGTTGGGAACGGACTTGCCGCCGTCGACGTTCCTGCTGTAGTTGTCGATGCTCACGACCACGGGAAAGCCGAAGCTCCCGTTCTTGACGCTCGACCTGCACTGGACGTTGTAGTTCTTGCCGGCGACGCTCGTCAAGACGTTCTTGTCGCGGCTGTTGTTTCACCACACATCGACGTGGCTCCCGATGTCTGACCGGCTCTTGGTCTTCTCGACGCTCGGGCGGCACGGGAAGGGGACGCCGTCTCACTCGCACAGGGTCTGGTCGTTCGACTGGTCTCGGCTCTGGCACAAGGGGCTGCACCTGATGCACTTATCGTGGCTCCGGCTCTGGTAGGACCTGTTGTAGTGGGTCTCGTGGGTCAGGAAGTTGCTGAAGTGGTCTCAACACCCGCCGCTCCTGCGGTTCGGGCCGGTCAAGGGGACCGTCCACCACGACTTGCCGTTTCACCTACGGAAGACGCCGCCGTCGTAGCACTTGCTCTTCACCTAGCACTGTCGGCGGGTGACGCACCTCTGGCCGCACTTCTAGTGGCACCACCGGCCGCTTGTGTTATAGCTCCTCTGGCTCGTGTGGCTCGTCTTCGCCTTGCAGTAGGCCTAATAGGGGGTGGTGTTGATGTTGCGGCGGTAGTTGTTCATGTTGGTGCTGTAGCGGGACGACCTCGACCTGCTCGGAGACCACGACTTATCGATGCACTGGGGGTAGACGTAGCGGCTGTTCCTCATGTGGTTGTAGAAGGACTTCAAGCCGTCGCCGATGCACAGGCCGACCCCGTCTCACAAGGTGTTCCCGTCTTCGCGGGACCACGACGTCATGGACTCTCACGGGGACCACCTGTCTCGGTGGACGGACGAATCGTGGTTCAAGTGGTAGATGTTGTTGTACAAGACGCGGCCGAAGGTGCTCCCGCCGTCTCTGTCGACGGTCCCGCTGTCGCCGCCTGGGGTGCACTGGCTTCACCTCCCGTGGTCGAAGGACTGGCCGTAGTAGTCGACCCCGCTCCTCACGCGGTACTTCCCGTTCATGCCGTAGATGTGGTTTCACTCGGCCATGCACTTGACCTAGTTCCTCTTTTGGTTCGACTGGACTGACTTAAGCAGCTGTTAGTTGGAGACCTAATGTTTTAAACACTTTCTAACTGACCATAAGAATTGATACAACGAGGAAAATGCGATACACCTATGCGACGAAATTACGGAAACATAGTACGATAACGAAGGGCATACCGAAAGTAAAAGAGGAGGAACATATTTAGGACCAACGACAGAGAAATACTCCTCAACACCGGGCAACAGTCCGTTGCACCGCACCACACGTGACACAAACGACTGCGTTGGGGGTGACCAACCCCGTAACGGTGGTGGACAGTCGAGGAAAGGCCCTGAAAGCGAAAGGGGGAGGGATAACGGTGCCGCCTTGAGTAGCGGCGGACGGAACGGGCGACGACCTGTCCCCGAGCCGACAACCCGTGACTGTTAAGGCACCACAACAGCCCCTTTAGTAGCAGGAAAGGAACCGACGAGCGGACACAACGGTGGACCTAAGACGCGCCCTGCAGGAAGACGATGCAGGGAAGCCGGGAGTTAGGTCGCCTGGAAGGAAGGGCGCCGGACGACGGCCGAGACGCCGGAGAAGGCGCAGAAGCGGAAGCGGGAGTCTGCTCAGCCTAGAGGGAAACCCGGCGGAGGGGCGGACCGACGAGCTCTCTAGCCCACCGTAGGGACACTGGGGAGGGGTCACGGAGAGGACCGGGACCTTCAACGGTGAGGTCACGGGTGGTCGGAACAGGATTATTTTAATTCAACGTAGTAAAACAGACTGATCCACAGGAAGATATTATAATACCCCACCTCCCCCCACCATACCTCGTTCCCCGTTCAACCCTTCTGTTGGACATCCCGGACGCCCCAGATAACCCTTGGTTCGACCTCACGTCACCGTGTTAGAACCGAGTGACGTTAGAGGCGGAGGACCCAAGTTCGCTAAGAGGACGGAGTCGGAGGGCTCAACAACCCTAAGGTCCGTACGTACTGGTCCGAGTCGATTAAAAACAAAAAAACCATCTCTGCCCCAAAGTGGTATAACCGGTCCGACCAGAGGTTGAGGATTAGAGTCCACTAGATGGGTGGAACCGGAGGGTTTAACGACCCTAATGTCCGCACTTGGTGACGAGGGAAGGGACAGGAATCCTTGGGGATCACTACCTCAACCGGTGAGGGAGAGACGCGCGAGCGAGCGAGTGACTCCGGCCCGCTGGTTTCCAGCGGGCTGCGGGCCCGAAACGGGCCCGCCGGAGTCACTCGCTCGCTCGCGCGTCGACGGACGTCCCCGCGGACTACGCCATAAAAGAGGAATGCGTAGACACGCCATAAAGTGTGGCGTATGCAGTTTCGTTGGTATCATGCGCGGGACATCGCCGCGTAATTCGCGCCGCCCACACCACCAATGCGCGTCGCACTGGCGATGTGAACGGTCGCGGAATCGCGGGCGAGGAAAGCGAAAGAAGGGAAGGAAAGAGCGGTGCAAGCGGCCGAAAGGGGCAGTTCGAGATTTAGCCCCCGAGGGAAATCCCAAGGCTAAATCACGAAATGCCGTGGAGCTGGGGTTTTTTGAACTAAACCCACTACCAAGTGCATCACCCGGTAGCGGGACTATCTGCCAAAAAGCGGGAAACTGCAACCTCAGGTGCAAGAAATTATCACCTGAGAACAAGGTTTGACCTTGTTGTGAGTTGAGATAGAGCCCGATAAGAAAACTAAATATTCCCTAAAACGGCTAAAGCCAGATAACCAATTTTTTACTCGACTAAATTGTTTTTAAATTGCGCTTAAAATTGTTTTATAATTGCAAATGTTAAAATACCACGTGAGAGTCATGTTAGACGAGACTACGGCGTATCAATTCGGTCGGGGCTGTGGGCGGTTGTGGGCGACTGCGCGGGACTGCCCGAACAGACGAGGGCCGTAGGCGAATGTCTGTTCGACACTGGCAGAGGCCCTCGACGTACACAGTCTCCAAAAGTGGCAGTAGTGGCTTTGCGCGCTCTGCTTTCCCGGAGCACTATGCGGATAAAAATATCCAATTACAGTACTATTATTACCAAAGAATCTGCAGTCCACCGTGAAAAGCCCCTTTACACGCGCCTTGGGGATAAACAAATAAAAAGATTTATGTAAGTTTATACATAGGCGAGTACTCTGTTATTGGGACTATTTACGAAGTTATTATAACTTTTTCCTTCTCATACTCATAAGTTGTAAAGGCACAGCGGGAATAAGGGAAAAAACGCCGTAAAACGGAAGGACAAAAACGAGTGGGTCTTTGCGACCACTTTCATTTTCTACGACTTCTAGTCAACCCACGTGCTCACCCAATGTAGCTTGACCTAGAGTTGTCGCCATTCTAGGAACTCTCAAAAGCGGGGCTTCTTGCAAAAGGTTACTACTCGTGAAAATTTCAAGACGATACACCGCGCCATAATAGGGCATAACTGCGGCCCGTTCTCGTTGAGCCAGCGGCGTATGTGATAAGAGTCTTACTGAACCAACTCATGAGTGGTCAGTGTCTTTTCGTAGAATGCCTACCGTACTGTCATTCTCTTAATACGTCACGACGGTATTGGTACTCACTATTGTGACGCCGGTTGAATGAAGACTGTTGCTAGCCTCCTGGCTTCCTCGATTGGCGAAAAAACGTGTTGTACCCCCTAGTACATTGAGCGGAACTAGCAACCCTTGGCCTCGACTTACTTCGGTATGGTTTGCTGCTCGCACTGTGGTGCTACGGACATCGTTACCGTTGTTGCAACGCGTTTGATAATTGACCGCTTGATGAATGAGATCGAAGGGCCGTTGTTAATTATCTGACCTACCTCCGCCTATTTCAACGTCCTGGTGAAGACGCGAGCCGGGAAGGCCGACCGACCAAATAACGACTATTTAGACCTCGGCCACTCGCACCCAGAGCGCCATAGTAACGTCGTGACCCCGGTCTACCATTCGGGAGGGCATAGCATCAATAGATGTGCTGCCCCTCAGTCCGTTGATACCTACTTGCTTTATCTGTCTAGCGACTCTATCCACGGAGTGACTAATTCGTAACCATTGACAGTCTGGTTCAAATGAGTATATATGAAATCTAACTAAATTTTGAAGTAAAAATTAAATTTTCCTAGATCCACTTCTAGGAAAAACTATTAGAGTACTGGTTTTAGGGAATTGCACTCAAAAGCAAGGTGACTCGCAGTCTGGGGCATCTTTTCTAGTTTCCTAGAAGAACTCTAGGAAAAAAAGACGCGCATTAGACGACGAACGTTTGTTTTTTTGGTGGCGATGGTCGCCACCAAACAAACGGCCTAGTTCTCGATGGTTGAGAAAAAGGCTTCCATTGACCGAAGTCGTCTCGCGTCTATGGTTTATGACAAGAAGATCACATCGGCATCAATCCGGTGGTGAAGTTCTTGAGACATCGTGGCGGATGTATGGAGCGAGACGATTAGGACAATGGTCACCGACGACGGTCACCGCTATTCAGCACAGAATGGCCCAACCTGAGTTCTGCTATCAATGGCCTATTCCGCGTCGCCAGCCCGACTTGCCCCCCAAGCACGTGTGTCGGGTCGAACCTCGCTTGCTGGATGTGGCTTGACTCTATGGATGTCGCACTCGATACTCTTTCGCGGTGCGAAGGGCTTCCCTCTTTCCGCCTGTCCATAGGCCATTCGCCGTCCCAGCCTTGTCCTCTCGCGTGCTCCCTCGAAGGTCCCCCTTTGCGGACCATAGAAATATCAGGACAGCCCAAAGCGGTGGAGACTGAACTCGCAGCTAAAAACACTACGAGCAGTCCCCCCGCCTCGGATACCTTTTTGCGGTCGTTGCGCCGGAAAAATGCCAAGGACCGGAAAACGACCGGAAAACGAGTGTACA"
    it "Shouldn't find any primers; Sequence: whole plasmid; pos: 0" $
      fmap symbol . _seq' <$> designPrimer True wholePlasmid2 0 `shouldBe` Left "Bio.Tools.Sequence.Primers.Optimization: all primers designed from given position don't bind to target."
    it "Shouldn't find any primers; Sequence: whole plasmid; pos: 7" $
      fmap symbol . _seq' <$> designPrimer True wholePlasmid2 7 `shouldBe` Left "Bio.Tools.Sequence.Primers.Optimization: all primers designed from given position don't bind to target."
    it "Sequence: whole plasmid; pos: 10" $
      fmap symbol . _seq' <$> designPrimer True wholePlasmid2 10 `shouldBe` Right "CGACGCGCGAGCGAGCG"
