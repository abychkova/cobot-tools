module Bio.Tools.Sequence.Primers.Constants
 ( annealingTemp
 , bindingRate
 , eTgtEFoldRel
 , maxPrimerGCContent
 , maxPrimerLength
 , maxPrimerMeltingTemp
 , minPrimerGCContent
 , minPrimerLength
 , minPrimerMeltingTemp
 , topN
 ) where

-- | Minimum recommended primer's length.
--
minPrimerLength :: Int
minPrimerLength = 17

-- | Maximum recommended primer's length.
--
maxPrimerLength :: Int
maxPrimerLength = 35

-- | Temperature under which annealing of primers happens.
--
annealingTemp :: Double
annealingTemp = 62

-- | Minimum recommended primer's melting temperature.
--
minPrimerMeltingTemp :: Int
minPrimerMeltingTemp = 55

-- | Maximum recommended primer's melting temperature.
--
maxPrimerMeltingTemp :: Int
maxPrimerMeltingTemp = 75

-- | Minimum recommended primer's GC-content.
--
minPrimerGCContent :: Float
minPrimerGCContent = 40

-- | Maximum recommended primer's GC-content.
--
maxPrimerGCContent :: Float
maxPrimerGCContent = 60

-- | Rate of nucleotides from primer that should bind to target area of source sequence.
--
bindingRate :: Float
bindingRate = 0.75

-- | Minimum allowed value of primer's tagret interaction energy to
-- primer's folding energy relation.
--
eTgtEFoldRel :: Float
eTgtEFoldRel = 2.7

-- | Number of candidates that are taken if none of candidates satisfies
-- GC-content condition.
--
topN :: Int
topN = 10
