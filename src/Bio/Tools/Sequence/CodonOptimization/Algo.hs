module Bio.Tools.Sequence.CodonOptimization.Algo
    ( optimizeAA
    , optimizeDNA
    , score
    , scoreSequence
    , scoreCmp
    ) where

import           Bio.NucleicAcid.Nucleotide                     (symbol)
import           Bio.NucleicAcid.Nucleotide.Type                (DNA (..))
import           Bio.Protein.AminoAcid.Type                     (AA (..))
import           Bio.Tools.Sequence.CodonOptimization.Constants (ak2Codon, ak2MaxFrequCodon,
                                                                 codon2ak,
                                                                 codonFrequencies, motiveScoreWindow, defaultMotiveScore, 
                                                                 forbiddenMotiveScore)
import           Bio.Tools.Sequence.CodonOptimization.Types     (CodonScoreConfig (..),
                                                                 standardTemperature)
import           Bio.Tools.Sequence.ViennaRNA.Fold              (fold)
import           Data.List                                      (foldl',
                                                                 maximumBy,
                                                                 take)
import           Data.Map                                       as Map (lookup)
import           Data.Maybe                                     (fromMaybe)
import           Text.Regex.TDFA                                ((=~))

-- | 'optimizeDNA' function does translation from [DNA] to [AA] and then calls 'optimizeAA'
optimizeDNA :: CodonScoreConfig -- ^ Config data object. Contains main parameters of codon-optimization and all parameters for scoring function
            -> [DNA]            -- ^ Initial, not optimized nucleotide sequence
            -> [DNA]            -- ^ Result, optimized nucleotide sequence
optimizeDNA cfg dna = optimizeAA cfg (translate dna)
  where
    translate :: [DNA] -> [AA]
    translate [] = []
    translate dnaSeq =
        case Map.lookup (take 3 dnaSeq) codon2ak of
            Just ak -> ak : translate (drop 3 dnaSeq)
            _       -> error $ "Unknown codon: " ++ show (take 3 dnaSeq)

-- | 'optimizeAA' function does codon-optimisation for incoming amino-acid sequence.
-- Incoming amino-acid sequence transformed to nucleotide sequence and optimized used the codon-optimization algorithm.
-- Algorithm described here doi: 10.1007/s11693-010-9062-3
optimizeAA :: CodonScoreConfig  -- ^ Config data object. Contains main parameters of codon-optimization and all parameters for scoring function
           -> [AA]              -- ^ Initial, not optimized amino-acid sequence
           -> [DNA]             -- ^ Result, optimized nucleotide sequence
optimizeAA cfg@(CodonScoreConfig organism initLen winLen _ _ _ _ _ _ _ _ _ _) aa = foldl' concatByScore initial variants
  where
    lenAA = length aa
    variants = generateVariants (drop initLen aa) winLen
    fequCodonsMap = ak2MaxFrequCodon organism
    initial = concatMap (\ak -> maybe "" fst (Map.lookup ak fequCodonsMap)) (take initLen aa)

    -- | 'concatByScore' function gets maximum by score variable string and then concat it to result string
    concatByScore :: [DNA]   -- ^ initial string
                  -> [[DNA]] -- ^ list of variable string
                  -> [DNA]   -- ^ result string
    concatByScore result vars
        | length result == 3 * (lenAA - winLen - 1) = result ++ maximumBy (scoreCmp cfg result) vars
        | otherwise = result ++ take 3 (maximumBy (scoreCmp cfg result) vars)


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

-- | 'scoreSequence' function calculates the average score for full sequence
scoreSequence :: CodonScoreConfig -> [DNA] -> Double
scoreSequence cnf@(CodonScoreConfig _ initLen winLen _ _ _ _ _ _ _ _ _ _) nkSequ = sum res / realToFrac (length res)
  where
    res = scr ((initLen + winLen + 1) * 3) []

    scr :: Int -> [Double] -> [Double]
    scr partLen acc | partLen > length nkSequ = acc
                    | otherwise = scr (partLen + winLen * 3) (score cnf (take partLen nkSequ) : acc)

-- | 'score' function gets scoring for incoming string.
-- Scoring function is a composite function of several scoring. More about scoring algorithm see here doi: 10.1007/s11693-010-9062-3
score :: CodonScoreConfig  -- ^ Config data object. Contains main parameters of codon-optimization and all parameters for scoring function
      -> [DNA]             -- ^ nucleotide sequence to score
      -> Double            -- ^ result score value
score (CodonScoreConfig organism _ winLen codonUsageWeight gcWeight gcFactor gcWindow rnaFoldingWeight
                        rnaFoldingFactor rnaFoldingWindow forbiddenDNAWeight gcContentDesired forbiddenRegexp) nkSequ =
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
    motiveScore sequ =
        if any (drop (length sequ - motiveScoreWindow) (map symbol sequ) =~) forbiddenRegexp
            then forbiddenMotiveScore
            else defaultMotiveScore


    -- | 'scoreRnaf' counts energy of RNA folding
    scoreRnaf :: [DNA] -> Int
    scoreRnaf sequ = truncate $ rnaFoldingWeight * (abs result ** rnaFoldingFactor)
      where
        result = fst $ fold standardTemperature (drop (length sequ - rnaFoldingWindow) sequ)

-- | 'scoreCmp' is compare function for two strings using 'score' function
scoreCmp :: CodonScoreConfig -> [DNA] -> [DNA] -> [DNA] -> Ordering
scoreCmp cfg optimized str1 str2 = compare score1 score2
  where
    score1 = score cfg (optimized ++ str1)
    score2 = score cfg (optimized ++ str2)
