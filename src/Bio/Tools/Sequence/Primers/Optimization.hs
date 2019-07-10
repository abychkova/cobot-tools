{-# LANGUAGE ViewPatterns #-}

module Bio.Tools.Sequence.Primers.Optimization
 ( designPrimer
 ) where

import           Bio.Chain                            (fromList)
import           Bio.Chain.Alignment                  (AffineGap (..),
                                                       AffineGap2,
                                                       AlignmentResult (..),
                                                       LocalAlignment (..),
                                                       Operation (..), align)
import           Bio.NucleicAcid.Chain                (NucleicAcidChain (..))
import           Bio.NucleicAcid.Nucleotide           (Complementary (..),
                                                       DNA (..), symbol)
import           Bio.Tools.Sequence.Primers.Constants as Constants (annealingTemp,
                                                                    bindingRate,
                                                                    eTgtEFoldRel,
                                                                    maxPrimerGCContent,
                                                                    maxPrimerLength,
                                                                    maxPrimerMeltingTemp,
                                                                    minPrimerGCContent,
                                                                    minPrimerLength,
                                                                    minPrimerMeltingTemp,
                                                                    topN)
import           Bio.Tools.Sequence.Primers.Types     (Primer,
                                                       ScoredPrimer (..), eTgt,
                                                       eTgtEFold, gcContent,
                                                       meltingTemp, seq')
import           Bio.Tools.Sequence.ViennaRNA.Cofold  (cofold)
import           Bio.Tools.Sequence.ViennaRNA.Fold    (fold)
import           Control.Lens                         ((&), (.~), (^.))
import           Control.Monad                        (when)
import           Control.Monad.Except                 (MonadError, throwError)
import           Data.List                            (sortOn)
import           Data.Maybe                           (catMaybes)
import           Data.Text                            (Text)

-- | Given 'DNA' sequence and position in that sequence designs forward primer
-- for that sequence. Primer will start at the given position. @isCyclic@ marks
-- whether the sequence that we design primer for is cyclic or not.
--
-- Flag @isCyclic@ also defines type of algorithm that will be used to design primer.
--
-- If @isCyclic@ is set to False, ViennaRNA will be used to check that primer
-- has no off-target interactions with given sequence. This is pretty accurate method,
-- but when it is used to process long sequences (more than 1000 bps), it's quite slow.
--
-- If @isCyclic@ is set to True, algorithm that uses several heuristics will be
-- used to check that primer has no off-target interactions with given sequence.
-- This method is less accurate than using ViennaRNA to calculate energy of primer's interaction
-- with sequence, but much more faster.
--
-- We did such segregation of algorithms, because cyclic DNA sequences (plasmids)
-- are very long and non-cyclic sequences that we work with are no longer than 1000 bps.
--
designPrimer :: MonadError Text m => Bool -> [DNA] -> Int -> m ScoredPrimer
designPrimer isCyclic dna' pos = do
    when (pos >= length dna') $ throwError "Bio.Tools.Sequence.Primers.Optimization: given position is out of range."

    withE <- energyFilter . gcContentFilter . lengthFilter $ candidates

    case gcOnEndFilter withE of
      [] -> throwError badPrimersError
      l  -> pure $ head $ sortOn (fmap negate . (^. gcContent)) l
  where
    dna        = if isCyclic then dna' <> take Constants.maxPrimerLength dna' else dna'
    candidates = genTemperatureCandidates (drop pos dna) []

    badPrimersError :: Text
    badPrimersError = "Bio.Tools.Sequence.Primers.Optimization: all primers designed from given position are seriously flawed."

    -- | Generates candidate primers with melting temperature in needed range.
    --
    genTemperatureCandidates :: [DNA] -> [DNA] -> [ScoredPrimer]
    genTemperatureCandidates [] _          = []
    genTemperatureCandidates (x : xs) base = res
      where
        newBase        = base <> [x]
        newPrimer      = toScored newBase
        primerWithTemp = calcMeltingTemp newPrimer

        Just mTemp = primerWithTemp ^. meltingTemp

        res | mTemp < Constants.minPrimerMeltingTemp = genTemperatureCandidates xs newBase
            | mTemp > Constants.maxPrimerMeltingTemp = []
            | otherwise                              = primerWithTemp : genTemperatureCandidates xs newBase

        toScored :: Primer -> ScoredPrimer
        toScored p = ScoredPrimer p Nothing Nothing Nothing Nothing

    -- | Filters 'ScoredPrimer's based on their length.
    --
    lengthFilter :: [ScoredPrimer] -> [ScoredPrimer]
    lengthFilter = filter lengthPred
      where
        lengthPred :: ScoredPrimer -> Bool
        lengthPred sp = Constants.minPrimerLength <= ls && ls <= Constants.maxPrimerLength
          where
            ls = length $ sp ^. seq'

    -- | Leaves only 'ScoredPrimer's that have GC on their 3' end.
    -- If no such primers are found, filtering doesn't happen.
    --
    gcOnEndFilter :: [ScoredPrimer] -> [ScoredPrimer]
    gcOnEndFilter sps | null filtered = sps
                      | otherwise     = filtered
      where
        filtered = filter gcOnEndPred sps

        gcOnEndPred :: ScoredPrimer -> Bool
        gcOnEndPred = (`elem` [DC, DG]) . last . (^. seq')

    -- | Leaves primers whose GC-content is in needed range. If no such primers
    -- are found, leaves @Constants.topN@ primers with GC-content nearest to needed range.
    --
    gcContentFilter :: [ScoredPrimer] -> [ScoredPrimer]
    gcContentFilter sps | null filtered = take Constants.topN sorted
                        | otherwise     = filtered
      where
        primersWithGC = fmap calcGCContent sps
        filtered      = filter gcContentPred primersWithGC

        sorted = sortOn scoringFunc primersWithGC

        gcContentPred :: ScoredPrimer -> Bool
        gcContentPred sp = Constants.minPrimerGCContent <= gcc && gcc <= Constants.maxPrimerGCContent
          where
            Just gcc = sp ^. gcContent

        scoringFunc :: ScoredPrimer -> Float
        scoringFunc sp = min (abs $ Constants.minPrimerGCContent - gcc) (abs $ Constants.maxPrimerGCContent - gcc)
          where
            Just gcc = sp ^. gcContent

    -- | Filters primers using energy characteristics.
    --
    energyFilter :: MonadError Text m => [ScoredPrimer] -> m [ScoredPrimer]
    energyFilter sps = fmap eTgtEFoldFilter . eTgtFilter $ sps

    -- | Leaves only primers whose relation of energy of interaction with target
    -- to energy of forming a secondary structure is higher then @Constants.eTgtEFoldRel@.
    --
    eTgtEFoldFilter :: [ScoredPrimer] -> [ScoredPrimer]
    eTgtEFoldFilter s = filter eTgtEFoldPred . fmap calcTargetFold $ s
      where
        eTgtEFoldPred :: ScoredPrimer -> Bool
        eTgtEFoldPred sp = etef >= Constants.eTgtEFoldRel
          where
            Just etef = sp ^. eTgtEFold

    -- | Leaves only primers that bind to target on source sequence.
    --
    eTgtFilter :: MonadError Text m => [ScoredPrimer] -> m [ScoredPrimer]
    eTgtFilter s | null res  = throwError badPositionError
                 | otherwise = pure res
      where
        res = catMaybes $ fmap (calcTarget isCyclic dna pos) s

        badPositionError :: Text
        badPositionError = "Bio.Tools.Sequence.Primers.Optimization: all primers designed from given position don't bind to target."

-- | Calculates melting temperature for 'ScoredPrimer'.
-- Temperature is calculated using the following formula: 4 * (C + G) + 2 * (A + T).
--
calcMeltingTemp :: ScoredPrimer -> ScoredPrimer
calcMeltingTemp sp = sp & meltingTemp .~ Just temperature
  where
    temperature = sum $ fmap tempForNuc (sp ^. seq')

    tempForNuc :: DNA -> Int
    tempForNuc DC = 4
    tempForNuc DG = 4
    tempForNuc DA = 2
    tempForNuc DT = 2

-- | Calculates GC-content for 'ScoredPrimer'.
-- GC-content is calulcated using the following formula: (G + C) / (A + T + G + C).
--
calcGCContent :: ScoredPrimer -> ScoredPrimer
calcGCContent sp = sp & gcContent .~ Just content
  where
    nGC     = length $ filter (\x -> x == DC || x == DG) $ sp ^. seq'
    content = fromIntegral nGC / fromIntegral (length $ sp ^. seq')

-- | Calculates relation of 'ScoredPrimers's target interaction energy to
-- its folding energy.
--
calcTargetFold :: ScoredPrimer -> ScoredPrimer
calcTargetFold sp = sp & eTgtEFold .~ Just tgtToFold
  where
    foldingEnergy = fst $ fold Constants.annealingTemp $ sp ^. seq'

    -- @calcTargetFold@ is used only after 'ScoredPrimer's target interaction energy
    -- has been calculated
    Just tgt  = sp ^. eTgt
    tgtToFold = abs $ tgt / toEps foldingEnergy

-- | Calculates energy of interaction of given 'DNA' and given 'DNA' sequence's
-- complementary strand.
--
cofoldEnergy :: [DNA] -> [DNA] -> Float
cofoldEnergy sp s = abs $ fst $ cofold Constants.annealingTemp (reverse sp, fmap cNA s)

-- | Local alignment algorithm that aligns arguments based on their complementarity
-- and doesn't allow gaps on the second argument of the alignment (that is primer in our terms).
--
-- We don't allow gaps on the primer during alignment, because this situation
-- is not physical: it means that hairpins can appear on the plasmid.
--
localAlignment :: LocalAlignment AffineGap2 DNA DNA
localAlignment = LocalAlignment complementScoring (plasmidGap, primerGap)
  where
    plasmidGap :: AffineGap
    plasmidGap = AffineGap (-5) (-1)

    primerGap :: AffineGap
    primerGap = AffineGap (-1000) (-1000)

    complementScoring :: DNA -> DNA -> Int
    complementScoring DA (cNA -> DT) = 3
    complementScoring DT (cNA -> DA) = 3
    complementScoring DC (cNA -> DG) = 5
    complementScoring DG (cNA -> DC) = 5
    complementScoring _ _            = -3

-- | Calculates energy of interaction of primer @sp@ with sequence @dna@.
-- If @sp@ doesn't bind to target area on @dna@ (target area starts at @sourceInd@),
-- then Nothing is returned.
--
-- There are different algorithms to check off-target interaction depending of whether
-- @dna@ is cyclic or not. It is defined by the first parameter.
--
calcTarget :: Bool -> [DNA] -> Int -> ScoredPrimer -> Maybe ScoredPrimer
-- algorithm for cyclic sequences
calcTarget True dna sourceInd sp | null offTargets || tgtE > maxOffTarget = res
                                 | otherwise                              = Nothing
  where
    primerSeq     = sp ^. seq'
    primerLen     = length primerSeq

    -- position in the @dna@, where primer reaches half of its length
    primerLenHalf = sourceInd + primerLen `div` 2

    -- cut all @matchingSeqs@ into k-mers of all sizes in range [@primerLenHalf@; @primerLen@].
    offTargets = concatMap seqToKMers matchingSeqs

    tgtE         = cofoldEnergy primerSeq primerSeq
    maxOffTarget = maximum $ fmap (cofoldEnergy primerSeq) offTargets

    res = pure $ sp & eTgt .~ Just tgtE

    seqToKMers :: [DNA] -> [[DNA]]
    seqToKMers s = concatMap (toKMers s) [primerLenHalf .. primerLen]
      where
        toKMers :: [a] -> Int -> [[a]]
        toKMers l k | length l < k = []
                    | otherwise    = take k l : toKMers (tail l) k

    -- | All parts of @dna@ that align the best on given @sp@ using algorithm @localAlignment@.
    -- We consider complementary strands of these parts to be most energetically preferrable
    -- off-targets for our primer.
    --
    matchingSeqs :: [[DNA]]
    matchingSeqs = fmap matchingSeq [invertedDna, dnaRevComp]
      where
        -- we invert the cyclic @dna@ in such way that target area is being cut in half,
        -- so that it won't be considered as off-target interaction
        invertedDna | [x]    <- invertPoints = drop x dna <> drop Constants.maxPrimerLength (take x dna)
                    | [x, y] <- invertPoints = take (y - x) . drop x $ dna
                    | otherwise              = error "This branch is never visited."

        -- here we check that @primerLenHalf@ is not in the part of sequence,
        -- that was appended at the end to create cyclic sequence.
        -- If it's not in this part, invert at @primerLenHalf@.
        -- Otherwise consider everything in between @modPos@ and (@modPos@ + (length @dna@ - @Constants.maxPrimerLength@))
        -- as inverted sequence.
        modPos       = primerLenHalf `mod` (length dna - Constants.maxPrimerLength)
        invertPoints | modPos >= Constants.maxPrimerLength = [primerLenHalf]
                     | otherwise                           = [modPos, modPos + (length dna - Constants.maxPrimerLength)]

        -- off-target interaction could happen in reversed direction
        dnaRevComp = reverse $ fmap cNA dna

        matchingSeq :: [DNA] -> [DNA]
        matchingSeq dna' = res'
          where
            ar = alignment $ alignmentFunc dna' primerSeq

            traceStart = toCoord $ last ar
            traceEnd   = toCoord $ head ar

            -- (traceStart, traceEnd) is inclusive range that describes off-target area
            -- in @dna'@. We want to extend that range by @primerLen@ in both directions
            -- to consider different possibilities
            (l, r)  = (max 0 (traceEnd - primerLen), min (length dna' - 1) (traceStart + primerLen))
            res'    = take (r - l + 1) $ drop l dna'

        toCoord :: Operation Int Int -> Int
        toCoord = getI

        alignmentFunc :: [DNA] -> [DNA] -> AlignmentResult (NucleicAcidChain Int DNA) (NucleicAcidChain Int DNA)
        alignmentFunc plasmid primer = alRes
          where
            plasmidC = toChain plasmid
            primerC  = toChain primer

            alRes = align localAlignment plasmidC primerC

            toChain :: [DNA] -> NucleicAcidChain Int DNA
            toChain = NucleicAcidChain . fromList
-- algorithm for linear sequences
calcTarget False dna sourceInd sp | Just e <- primerTargetEnergy = pure $ sp & eTgt .~ Just e
                                  | otherwise = Nothing
  where
    primerTargetEnergy :: Maybe Float
    primerTargetEnergy | abs e < revE = Nothing
                       | Just cnt <- actualBindingRateM, checkBindingRate cnt = Just e
                       | otherwise    = Nothing
      where
        primerSeq = sp ^. seq'
        spStr     = symbol <$> primerSeq

        compPrimerSeq = fmap cNA dna

        (e, bindStr) = cofold Constants.annealingTemp (reverse primerSeq, compPrimerSeq)
        primerLength = length spStr

        -- energy of primer's interaction revesed dna strand
        revE = cofoldEnergy primerSeq (reverse $ fmap cNA dna)

        -- inclusive range in which target binding site is being contained in @bindStr@
        (lInd, rInd) = (primerLength + sourceInd, lInd + primerLength - 1)

        bindWithInds       = zip bindStr [0..]
        actualBindingRateM = calcMatchesInBindingSite [] 0 bindWithInds

        -- | Calculates number of nucleotides in target binding site that bind
        -- to primer. Also checks that last nucleotide of primer binds to last nucleotide
        -- of target binding site. If this condition is not satisfied, Nothing is returned.
        --
        -- @bindStr@ is used for these calculations. @bindStr@ represents interaction between primer
        -- and sequence in form of a dot plot. Dot plot also contains balanced bracket sequence
        -- that shows how nucleotides bind to each other.
        --
        calcMatchesInBindingSite :: [Int] -> Int -> [(Char, Int)] -> Maybe Int
        calcMatchesInBindingSite _ cnt []           = Just cnt
        calcMatchesInBindingSite stack cnt (x : xs) -- if we encounter close bracket, then we check that its position is in target binding site
                                                    -- and position of open bracket that corresponds to it is in primer
                                                    | c == closeBracket, lInd <= i && i < rInd && i' < primerLength = calcMatchesInBindingSite l (cnt + 1) xs
                                                    -- close bracket in target corresponds to non-primer interaction, we skip this bracket
                                                    | c == closeBracket, lInd <= i && i < rInd = calcMatchesInBindingSite l cnt xs
                                                    -- next two conditions check that 3' end of primer binds to end of the target binding site
                                                    | c == closeBracket, i == rInd, i' == 0 = Just $ cnt + 1
                                                    | c == closeBracket, i == rInd, i' /= 0 = Nothing
                                                    -- close bracket is out of target range, we skip it
                                                    | c == closeBracket = calcMatchesInBindingSite l cnt xs
                                                    -- open brackets are put on stack
                                                    | c == openBracket = calcMatchesInBindingSite (i : stack) cnt xs
                                                    -- we ignore dots
                                                    | otherwise        = calcMatchesInBindingSite stack cnt xs
          where
            (c, i) = x

            -- since the bracket sequence is balanced, we will get to this pattern-matching only if
            -- there is something on top of the stack
            (i' : l) = stack

        openBracket :: Char
        openBracket = '('

        closeBracket :: Char
        closeBracket = ')'

        -- | Checks that not less then @Constants.bindingRate@ * 100 precents of primer's nucleotides
        -- interact with target.
        --
        checkBindingRate :: Int -> Bool
        checkBindingRate cnt = fromIntegral cnt / fromIntegral primerLength >= Constants.bindingRate

--------------------------------------------------------------------------------
-- Utility functions.
--------------------------------------------------------------------------------

-- | This function is used to avoid dividing by zero.
--
toEps :: Float -> Float
toEps x | abs x < eps = eps
        | otherwise   = x
  where
    eps = 0.01
