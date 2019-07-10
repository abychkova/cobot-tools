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
