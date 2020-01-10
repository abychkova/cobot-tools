import           SpecPrimers
import           SpecCodonOptimization
import           OligoDesigner.SpecOligoDesignerSplitter
import           OligoDesigner.SpecOligoDesignerScorer
import           OligoDesigner.SpecUtils
import           SpecViennaRNA
import           System.IO
import           Test.Hspec

main :: IO ()
main = do
    hSetBuffering stdout NoBuffering
    hspec $ do
         -- Primers
--         testPrimers
--
--         -- ViennaRNA
--         foldTest
--         cofoldTest
--
--         -- CodonOptimization
--         codonOptimizationSpec

         -- OligoDesigner
         oligoDesignerSplitterSpec
         oligoDesignerScoreSpec
         utilsSpec

