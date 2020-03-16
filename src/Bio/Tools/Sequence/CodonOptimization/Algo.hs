module Bio.Tools.Sequence.CodonOptimization.Algo
    ( optimizeCodonForAA
    , optimizeCodonForDNA
    , scoreByWindow
    , score
    ) where

import           Bio.NucleicAcid.Nucleotide                     (symbol)
import           Bio.NucleicAcid.Nucleotide.Type                (DNA (..))
import           Bio.Protein.AminoAcid.Type                     (AA (..))
import           Bio.Tools.Sequence.CodonOptimization.Constants (ak2Codon, ak2MaxFrequCodon,
                                                                 codon2ak,
                                                                 codonFrequencies,
                                                                 defaultMotiveScore,
                                                                 forbiddenMotiveScore,
                                                                 motiveScoreWindow)
import           Bio.Tools.Sequence.CodonOptimization.Types     (CodonOptimizationConfig (..),
                                                                 standardTemperature)
import           Bio.Tools.Sequence.ViennaRNA.Fold              (fold)
import           Data.List                                      (foldl',
                                                                 maximumBy,
                                                                 take)
import           Data.Map                                       as Map (lookup)
import           Data.Maybe                                     (fromMaybe)
import           Text.Regex.TDFA                                (makeRegex, Regex, match)
import Debug.Trace

-- | optimizeCodonForDNA function does translation from [DNA] to [AA] and then calls 'optimizeCodonForAA'
optimizeCodonForDNA :: CodonOptimizationConfig -- ^ Config data object. Contains main parameters of codon-optimization and all parameters for scoring function
            -> [DNA]            -- ^ Initial, not optimized nucleotide sequence
            -> [DNA]            -- ^ Result, optimized nucleotide sequence
optimizeCodonForDNA cfg dna = optimizeCodonForAA cfg (translate dna)
  where
    translate :: [DNA] -> [AA]
    translate [] = []
    translate dnaSeq =
        case Map.lookup (take 3 dnaSeq) codon2ak of
            Just ak -> ak : translate (drop 3 dnaSeq)
            _       -> error $ "Unknown codon: " ++ show (take 3 dnaSeq)

-- | 'optimizeCodonForAA' function does codon-optimisation for incoming amino-acid sequence.
-- Incoming amino-acid sequence transformed to nucleotide sequence and optimized used the codon-optimization algorithm.
-- Algorithm described here doi: 10.1007/s11693-010-9062-3
optimizeCodonForAA :: CodonOptimizationConfig  -- ^ Config data object. Contains main parameters of codon-optimization and all parameters for scoring function
           -> [AA]              -- ^ Initial, not optimized amino-acid sequence
           -> [DNA]             -- ^ Result, optimized nucleotide sequence
optimizeCodonForAA cfg@(CodonOptimizationConfig organism initLen winLen _ _ _ _ _ _ _ _ _ forbiddenRegexp) aa =
    traceMarker "optimizeCodonForAA: line 47" $ foldl' concatByScore initial variants
  where
    regex = map makeRegex forbiddenRegexp :: [Regex]
    lenAA = length aa
    variants = generateVariants (drop initLen aa) winLen
    fequCodonsMap = ak2MaxFrequCodon organism
    initial = concatMap (\ak -> maybe "" fst (Map.lookup ak fequCodonsMap)) (take initLen aa)

    -- | 'concatByScore' function gets maximum by score variable string and then concat it to result string
    concatByScore :: [DNA]   -- ^ initial string
                  -> [[DNA]] -- ^ list of variable string
                  -> [DNA]   -- ^ result string
    concatByScore result vars
        | length result == 3 * (lenAA - winLen - 1) =  result ++ fst (maximumBy compareBySecond (map var2score vars))
        | otherwise = result ++ take 3 (fst (maximumBy compareBySecond (map var2score vars)))
      where
        var2score :: [DNA] -> ([DNA], Double)
        var2score dna = (dna, scoreByWindow cfg regex (result ++ dna))

        compareBySecond :: (a, Double) -> (a, Double) -> Ordering
        compareBySecond p1 p2 = compare (snd p1) (snd p2)


-- | 'generateVariants' function generates list of all possible variants of nucleotide sequence for amino-acid sequence.
-- It is just recursive execution of 'windowVariants' for all amino-acid sequence.
-- Example: generateVariants [PHE,TRP,GLU,MET] 2 => [[[DT,DT,DT,DT,DG,DG,DG,DA,DA], [DT,DT,DC,DT,DG,DG,DG,DA,DA], [DT,DT,DT,DT,DG,DG,DG,DA,DG], [DT,DT,DC,DT,DG,DG,DG,DA,DG]],
--                                        [[DT,DG,DG,DG,DA,DA,DA,DT,DG], [DT,DG,DG,DG,DA,DG,DA,DT,DG]]]
-- Returns empty list in case of empty incoming string
generateVariants :: [AA]      -- ^ amino-acid sequence
                 -> Int       -- ^ length of window. means how much amino-acid from the right side side will be taken during scoring variant for one codon
                 -> [[[DNA]]] -- ^ result list. for each amino-acid now there is list of all variants for nucleotide sequence
generateVariants [] _ = []
generateVariants aa winLen
    | length aa == winLen = []
    | otherwise = windowVariants aa winLen : generateVariants (drop 1 aa) winLen

-- | 'windowVariants' function generates list of all possible variants of nucleotide sequence for amino-acid window.
-- Example: windowVariants [PHE,TRP,GLU,MET] 2 => [[DT,DT,DT,DT,DG,DG,DG,DA,DA], [DT,DT,DC,DT,DG,DG,DG,DA,DA], [DT,DT,DT,DT,DG,DG,DG,DA,DG], [DT,DT,DC,DT,DG,DG,DG,DA,DG]]
-- Returns empty list in case of empty incoming string
windowVariants :: [AA] -> Int -> [[DNA]]
windowVariants sequ winLen = map concat . mapM getCodons . take (winLen + 1) $ sequ

-- | 'getCodons' function gets list of codons for amino-acid
-- Example: 'getCodons' PRO => [[DC,DC,DT], [DC,DC,DC], [DC,DC,DA], [DC,DC,DG]]
-- Returns empty list in case of unknown amino-acid
getCodons :: AA -> [[DNA]]
getCodons ak = fromMaybe [] (Map.lookup ak ak2Codon)

-- | 'score' function calculates the average score for full sequence
score :: CodonOptimizationConfig -> [DNA] -> Double
score cnf@(CodonOptimizationConfig _ initLen winLen _ _ _ _ _ _ _ _ _ forbiddenRegexp) nkSequ =
    sum res / realToFrac (length res)
  where
    regex = map makeRegex forbiddenRegexp :: [Regex]
    res = scr ((initLen + winLen + 1) * 3) []

    scr :: Int -> [Double] -> [Double]
    scr partLen acc | partLen > length nkSequ = acc
                    | otherwise = scr (partLen + winLen * 3) (scoreByWindow cnf regex (take partLen nkSequ) : acc)

-- | 'scoreByWindow' function gets scoring for incoming string.
-- Scoring function is a composite function of several scoring. More about scoring algorithm see here doi: 10.1007/s11693-010-9062-3
scoreByWindow :: CodonOptimizationConfig  -- ^ Config data object. Contains main parameters of codon-optimization and all parameters for scoring function
      -> [Regex]           -- ^ compiled regexes for forbidden sequences
      -> [DNA]             -- ^ nucleotide sequence to score
      -> Double            -- ^ result score value
scoreByWindow (CodonOptimizationConfig organism _ winLen codonUsageWeight gcWeight gcFactor gcWindow rnaFoldingWeight
                        rnaFoldingFactor rnaFoldingWindow forbiddenDNAWeight gcContentDesired _) forbiddenRegexes nkSequ =
    scoreCU + scoreGC - scoreMT - realToFrac scoreRNAFold
  where
    sequLen = length nkSequ
    optimizedLen = sequLen - (winLen + 1) * 3
    scoreGC =
        if optimizedLen < gcWindow - optimizedLen  -- check if we have enough sequence for gcWindow
            then 0
            else realToFrac gcWeight * gcScore (drop (sequLen - gcWindow) nkSequ)
    scoreCU = realToFrac codonUsageWeight * codonUsage (drop optimizedLen nkSequ)
    scoreMT = realToFrac forbiddenDNAWeight * motiveScore nkSequ
    scoreRNAFold = scoreRnaf nkSequ

    -- | 'gcScore' function for the GC content.
    -- It is negatively counted absolute difference between the desired GC content and the GC content of a test sequence.
    gcScore :: [DNA] -> Double
    gcScore sequ = -abs (gc / (at + gc) * 100 - realToFrac gcContentDesired) ** gcFactor
      where
        gc = realToFrac $ length $ filter (\s -> s == DC || s == DG) sequ
        at = realToFrac $ length $ filter (\s -> s == DA || s == DT) sequ

    -- | 'codonUsage' function gets higest score for most frequently used codons
    codonUsage :: [DNA] -> Double
    codonUsage sequ = (cai ** (1 / codonCount)) * 100
      where
        codonCount = realToFrac (length sequ) / 3
        cai = countWeight sequ 1.0

        -- | 'countWeight' function is recursive counting weight for incoming string according to codon usage frequencies
        countWeight :: [DNA] -> Double -> Double
        countWeight [] acc = acc
        countWeight windowSeq acc = countWeight (drop 3 windowSeq) (acc * fromMaybe 0 (countWeightMb windowSeq))

        countWeightMb :: [DNA] -> Maybe Double
        countWeightMb str = do
            let codon = take 3 str
            codonFreq         <- Map.lookup codon (codonFrequencies organism)
            ak                <- Map.lookup codon codon2ak
            (_, codonMaxFreq) <- Map.lookup ak (ak2MaxFrequCodon organism)
            return $ codonFreq / codonMaxFreq

    -- | 'motiveScore' counts score for the occurrence of desired and unwanted DNA motifs.
    motiveScore :: [DNA] -> Double
    motiveScore sequ = res
      where
        dnaStr = map symbol (drop (length sequ - motiveScoreWindow) sequ)
        res = if any (`match` dnaStr) forbiddenRegexes
            then forbiddenMotiveScore
            else defaultMotiveScore

    -- | 'scoreRnaf' counts energy of RNA folding
    scoreRnaf :: [DNA] -> Int
    scoreRnaf sequ = truncate $ rnaFoldingWeight * (abs result ** rnaFoldingFactor)
      where
        result = fst $ fold standardTemperature (drop (length sequ - rnaFoldingWindow) sequ)
