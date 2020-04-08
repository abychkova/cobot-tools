module Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer(
    fixForbidden
   ,filterForbidden
) where

import Bio.NucleicAcid.Nucleotide (DNA)

import System.Random (StdGen)

import Control.Monad.State (State, evalState)
import Control.Monad.Except (Except, throwError)
import Text.Regex.TDFA (Regex, makeRegex, match, getAllMatches)
import Bio.Tools.Sequence.OligoDesigner.Utils.Prettifier (prettyDNA)
import Bio.Tools.Sequence.OligoDesigner.Utils.MutationUtils (mutate)
import Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils (getAAIndex)
import Bio.Tools.Sequence.OligoDesigner.Types (OligsDesignerInnerConfig(..))
import Bio.Tools.Sequence.CodonOptimization (forbiddenSequence, organism)
import Bio.Tools.Sequence.CodonOptimization.Types (Organism)
import Debug.Trace (trace)

fixForbidden :: StdGen -> OligsDesignerInnerConfig -> [DNA] -> Except String [DNA]
fixForbidden gen (OligsDesignerInnerConfig organism _ regexes _ maxIteration) dna = do
    let res = evalState (fixIterative organism regexes (maxIteration + 1) dna) gen
    if null res then throwError "cannot fix this shit" else return res

fixIterative :: Organism -> [Regex] -> Int -> [DNA] -> State StdGen [DNA]
fixIterative organism regexes 0 dna         = return []
fixIterative organism regexes iteration dna =
    case getPositions (prettyDNA dna) of
        []        -> return dna
        positions -> trace ("dna:" ++ prettyDNA dna) $ trace ("fix positions:" ++ show positions) $ fixPositions positions [dna]
  where
    fixPositions :: [(Int, Int)] -> [[DNA]] -> State StdGen [DNA]
    fixPositions [] results               = do
        let filtered = filterForbidden regexes results
        if null filtered then return [] else return $ head filtered
    fixPositions (position : xs) results = do
      variants <- sequence [mutate organism result position | result <- results]
      fixPositions xs (concat variants)

    getPositions :: String -> [(Int, Int)]
    getPositions dna = res
      where
        matches = concat [getAllMatches (regex `match` dna) :: [(Int, Int)] | regex <- regexes]
        res = [(getAAIndex begin, getAAIndex (begin + len)) | (begin, len) <- matches]

filterForbidden :: [Regex] -> [[DNA]] -> [[DNA]]
filterForbidden regexes = filter notMatch
  where
    notMatch :: [DNA] -> Bool
    notMatch dna = True `notElem` [regex `match` prettyDNA dna :: Bool | regex <- regexes]