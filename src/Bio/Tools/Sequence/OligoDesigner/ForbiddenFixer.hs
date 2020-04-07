module Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer(
    fixForbidden
) where

import Bio.NucleicAcid.Nucleotide (DNA)

import System.Random (StdGen)

import Control.Monad.State (State, evalState)
import Control.Monad.Except (Except, throwError)
import Text.Regex.TDFA (Regex, makeRegex, match, getAllMatches)
import Bio.Tools.Sequence.OligoDesigner.Prettifier (prettyDNA)
import Bio.Tools.Sequence.OligoDesigner.Utils (mutate, getAAIndex, notMatch)
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
fixIterative organismType regexes 0 dna         = return []
fixIterative organismType regexes iteration dna =
    case getPositions (prettyDNA dna) of
        []        -> return dna
        positions -> trace ("dna:" ++ prettyDNA dna) $ trace ("fix positions:" ++ show positions) $ fixPositions positions [dna]
  where
    fixPositions :: [(Int, Int)] -> [[DNA]] -> State StdGen [DNA]
    fixPositions [] results               = do
        let filtered = filter (notMatch regexes) results
        if null filtered then return [] else return $ head filtered
    fixPositions (position : xs) results = do
      variants <- sequence [mutate organismType [] result position | result <- results]
      fixPositions xs (concat variants)

    getPositions :: String -> [(Int, Int)]
    getPositions dna = res
      where
        matches = concat [getAllMatches (regex `match` dna) :: [(Int, Int)] | regex <- regexes]
        res = [(getAAIndex begin, getAAIndex (begin + len)) | (begin, len) <- matches]