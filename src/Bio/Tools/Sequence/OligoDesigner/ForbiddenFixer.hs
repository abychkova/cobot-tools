module Bio.Tools.Sequence.OligoDesigner.ForbiddenFixer
    ( fixForbidden
    , filterForbidden
    ) where

import Control.Monad.Except           (Except, throwError)
import Control.Monad.Trans            (lift)
import Control.Monad.Trans.State.Lazy (StateT)
import System.Random                  (StdGen)
import Text.Regex.TDFA                (Regex, getAllMatches, match)

import Bio.NucleicAcid.Nucleotide (DNA)

import Bio.Tools.Sequence.CodonOptimization.Types           (Organism)
import Bio.Tools.Sequence.OligoDesigner.Types               (OligsDesignerInnerConfig (..))
import Bio.Tools.Sequence.OligoDesigner.Utils.CommonUtils   (getAAIndex)
import Bio.Tools.Sequence.OligoDesigner.Utils.MutationUtils (mutate)
import Bio.Tools.Sequence.OligoDesigner.Utils.Prettifier    (prettyDNA)


fixForbidden :: OligsDesignerInnerConfig -> [DNA] -> StateT StdGen (Except String) [DNA]
fixForbidden (OligsDesignerInnerConfig organism _ regexes _ maxIteration) dna = do
    res <- fixIterative organism regexes (maxIteration + 1) dna
    if null res then throwError "cannot fix this shit" else return res

fixIterative :: Organism -> [Regex] -> Int -> [DNA] -> StateT StdGen (Except String) [DNA]
fixIterative _ _ 0 _  = return []
fixIterative organism regexes iteration dna = do
    forbiddenPositions <- lift $ getPositions (prettyDNA dna)
    case forbiddenPositions of
        []        -> return dna
        positions -> fixPositions positions [dna]
  where
    fixPositions :: [(Int, Int)] -> [[DNA]] -> StateT StdGen (Except String) [DNA]
    fixPositions [] results = do
        let filtered = filterForbidden regexes results
        if null filtered then fixIterative organism regexes (iteration - 1) dna else return $ head filtered
    fixPositions (position : xs) results = do
      variants <- traverse (\result -> mutate organism result position) results
      fixPositions xs (concat variants)

    getPositions :: String -> Except String [(Int, Int)]
    getPositions sequ = do
        let matches = concat [getAllMatches (regex `match` sequ) :: [(Int, Int)] | regex <- regexes]
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
    notMatch dna = not $ any (`match` prettyDNA dna) regexes
