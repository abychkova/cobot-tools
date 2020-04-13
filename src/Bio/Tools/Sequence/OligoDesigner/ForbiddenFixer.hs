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
import Control.Monad.Trans.State.Lazy (StateT, get)
import Control.Monad.Trans (lift)

fixForbidden :: OligsDesignerInnerConfig -> [DNA] -> StateT StdGen (Except String) [DNA]
fixForbidden (OligsDesignerInnerConfig organism _ regexes _ maxIteration) dna = do
    res <- fixIterative organism regexes (maxIteration + 1) dna
    if null res then throwError "cannot fix this shit" else return res

fixIterative :: Organism -> [Regex] -> Int -> [DNA] -> StateT StdGen (Except String) [DNA]
fixIterative organism regexes 0 dna         = return []
fixIterative organism regexes iteration dna = do
    positions <- lift $ getPositions (prettyDNA dna)
    case positions of
        []        -> return dna
        positions -> trace ("dna:" ++ prettyDNA dna) $ trace ("fix positions:" ++ show positions) $ fixPositions positions [dna]
  where
    fixPositions :: [(Int, Int)] -> [[DNA]] -> StateT StdGen (Except String) [DNA]
    fixPositions [] results = do
        let filtered = filterForbidden regexes results
        if null filtered then return [] else return $ head filtered
    fixPositions (position : xs) results = do
      variants <- sequence [mutate organism dna position | dna <- results]
      fixPositions xs (concat variants)

    getPositions :: String -> Except String [(Int, Int)]
    getPositions dna = do
        let matches = concat [getAllMatches (regex `match` dna) :: [(Int, Int)] | regex <- regexes]
        mapM toAAIndexes matches
      where
        toAAIndexes :: (Int, Int) -> Except String (Int, Int)
        toAAIndexes (begin, len) = do
            start <- getAAIndex begin
            finish <- getAAIndex (begin + len)
            return (start, finish)

filterForbidden :: [Regex] -> [[DNA]] -> [[DNA]]
filterForbidden regexes = filter notMatch
  where
    notMatch :: [DNA] -> Bool
    notMatch dna = True `notElem` [regex `match` prettyDNA dna :: Bool | regex <- regexes]