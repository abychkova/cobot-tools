import           SpecPrimers
import           SpecCodonOptimization
import           SpecOligoDesigner
import           SpecOligoDesignerScorer
import           SpecViennaRNA
import           System.IO
import           Test.Hspec

main :: IO ()
main = do
    hSetBuffering stdout NoBuffering
    hspec $ do
--         -- Primers
--         testPrimers
--
--         -- ViennaRNA
--         foldTest
--         cofoldTest
--
--         -- CodonOptimization
--         codonOptimizationSpec

         --OligoDesigner
--         oligoDesignerSpec
        oligoDesignerScoreSpec

