import           OligoDesigner.SpecScorer
import           OligoDesigner.SpecSplitter
import           OligoDesigner.SpecRNACofoldOptimizer
import           OligoDesigner.SpecGCContentOptimizer
import           OligoDesigner.SpecUtils
import           SpecCodonOptimization
import           SpecPrimers
import           SpecViennaRNA
import           System.IO
import           Test.Hspec

main :: IO ()
main = do
    hSetBuffering stdout NoBuffering
    hspec $ do
         -- Primers
         testPrimers

         -- ViennaRNA
         foldTest
         cofoldTest

         -- CodonOptimization
         codonOptimizationSpec

         -- OligoDesigner
         oligoDesignerSplitterSpec
         oligoDesignerScoreSpec
         utilsSpec
         optimizerSpec
         gcContentOptimizerSpec

