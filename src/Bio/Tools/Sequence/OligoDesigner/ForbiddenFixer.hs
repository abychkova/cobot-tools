module Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer(
    fixForbidden
) where

import Bio.NucleicAcid.Nucleotide (DNA)

import System.Random (StdGen)

import Control.Monad.State (State, evalState)
import Control.Monad.Except (Except, throwError)
import Text.Regex.TDFA (Regex, makeRegex, match)
import Bio.Tools.Sequence.OligoDesigner.Prettifier (prettyDNA)
import Bio.Tools.Sequence.OligoDesigner.Utils (mutate, getAANumber, notMatch)
import Bio.Tools.Sequence.OligoDesigner.Types (OligsDesignerConfig(..))
import Bio.Tools.Sequence.CodonOptimization (forbiddenSequence, organism)
import Bio.Tools.Sequence.CodonOptimization.Types (Organism)
import Debug.Trace (trace)

fixForbidden :: StdGen -> OligsDesignerConfig -> [Regex] -> [DNA] -> Except String [DNA]
fixForbidden gen conf regexes dna = do
    let organismType = organism $ codonOptimizationConfig conf
    let maxIteration = maxFixForbiddenIteration conf
    let res = evalState (fixIterative organismType regexes (maxIteration + 1) dna) gen
    if null res then throwError "cannot fix this shit" else return res

fixIterative :: Organism -> [Regex] -> Int -> [DNA] -> State StdGen [DNA]
fixIterative organismType regexes 0 dna = trace "zero iteration" $ return []
fixIterative organismType regexes iteration dna =
    case getPositions regexes (prettyDNA dna) of
        []        -> trace "empty positions" $ return dna
        positions -> fixPositions positions
  where
    fixPositions :: [(Int, Int)] -> State StdGen [DNA]
    fixPositions positions = do
        variants <- concat <$> sequence [mutate organismType regexes dna position | position <- positions]
        if null variants
            then fixIterative organismType regexes (iteration - 1) dna --рассчитываем, что в следующей мутации будут другие варианты
            else return $ head variants

getPositions :: [Regex] -> String -> [(Int, Int)]
getPositions regexes dna = res
  where
    matches = [regex `match` dna :: (String, String, String) | regex <- regexes]
    founded = filter (\(begin, found, final) -> not (null found)) matches
    res = [(getAANumber $ length begin, getAANumber (length begin + length final)) | (begin, found, final) <- founded]