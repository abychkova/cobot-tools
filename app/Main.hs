{-# LANGUAGE OverloadedStrings #-}

module Main where
import Bio.Tools.Sequence.OligoDesigner.Algo (designOligsAA, getRandomSeed)
import Control.Monad.Except (runExcept)
import Bio.Tools.Sequence.OligoDesigner.Types (OligsDesignerConfig(..), OligSet(..), Olig(..), OligSplitting(..))
import Bio.Tools.Sequence.OligoDesigner.Utils.Prettifier (prettyOligSet, prettyDNA)
import Data.Default (def)
import Debug.Trace (trace)
import           Control.Monad.IO.Class      (MonadIO)
import Bio.Tools.Sequence.OligoDesigner.Scorer (rnaScore, oligsGCContentDifference, gcContentScoreByOligs, commonScore)
import Data.Text (Text)
import System.Environment
import Data.String (fromString)
import           Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils       (assemble)
import Bio.NucleicAcid.Nucleotide.Type (DNA (..))
import Control.Monad.State (State, get, put, evalState)
import System.Random (StdGen, randomR)
import Bio.Protein.AminoAcid (AA)
import Bio.Tools.Sequence.CodonOptimization (CodonOptimizationConfig(..))
import Bio.Tools.Sequence.CodonOptimization.Types (Organism(..), defaultForbiddenRegexp)
import Data.Char (toUpper)
import Text.Regex.TDFA (Regex, makeRegex)
import Bio.Tools.Sequence.CodonOptimization.Algo (scoreByWindow)
import Bio.Tools.Sequence.OligoDesigner.Utils.Printer (printResult)

--for performance testing
main :: IO ()
main = do
    args <- getArgs
    seed <- getRandomSeed
    let dna = fromString $ map toUpper (head args)
--    let rnaFactor = read $ args !! 1
--    let oligsGCContentFactor = read $ args !! 2
--    let gcContentFactor = read $ args !! 3

    let coeffs = evalState randomCoeff seed
    let coeffs' = map (\((a, b, c), gen) -> ((realToFrac a / 10, realToFrac b / 10, realToFrac c / 10), gen)) coeffs
    runOpt coeffs' dna

  where
    runOpt :: [((Double, Double, Double), StdGen)] -> [AA] -> IO ()
    runOpt [] _ = return ()
    runOpt (coeffs : xs) dna = runOneOpt coeffs dna >> runOpt xs dna

    runOneOpt :: ((Double, Double, Double), StdGen) -> [AA] -> IO ()
    runOneOpt ((rnaFactor, oligsGCContentFactor, gcContentFactor), seed) dna = do
        let codonConf = CodonOptimizationConfig CHO 3 1 1 0.5 1.4 40 0.001 2.6 100 1 51 defaultForbiddenRegexp
        let conf = OligsDesignerConfig codonConf def 150 30
        let result = runExcept $ designOligsAA seed conf dna
        case result of
            Left err -> print ("error:" ++ err)
            Right value -> print ("oligs:" ++ prettyOligSet value) >> printResult conf "O" value

randomCoeff :: State StdGen [((Integer, Integer, Integer), StdGen)]
randomCoeff = tt 3 []
  where
    tt :: Integer -> [((Integer, Integer, Integer), StdGen)] -> State StdGen [((Integer, Integer, Integer), StdGen)]
    tt 0 acc = return acc
    tt count acc = do
        a <- random
        b <- random
        c <- random
        gen <- get
        if a + b + c - 10 == 0 then tt (count - 1) (((a, b, c), gen) : acc) else tt count acc

random :: State StdGen Integer
random = do
    gen <- get
    let (random, newGen) = randomR (0, 10) gen
    put newGen
    return random