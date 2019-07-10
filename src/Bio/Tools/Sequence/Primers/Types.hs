{-# LANGUAGE TemplateHaskell #-}

module Bio.Tools.Sequence.Primers.Types where

import           Bio.NucleicAcid.Nucleotide (DNA)
import           Control.Lens               (makeLenses)


-- | Primer is just an alias for sequence of nucleotides.
--
type Primer = [DNA]

-- | Primer with its characteristics.
--
data ScoredPrimer = ScoredPrimer { _seq'        :: Primer      -- ^ primer's sequence
                                 , _meltingTemp :: Maybe Int   -- ^ melting temperature
                                 , _gcContent   :: Maybe Float -- ^ GC-content
                                 , _eTgtEFold   :: Maybe Float -- ^ relation of energy of interaction with
                                                               --   target to energy of primer's forming
                                                               --   secondary structure
                                 , _eTgt        :: Maybe Float -- ^ energy of interaction with target if
                                                               --   primer binds to target. Otherwise set to 0
                                 }
  deriving (Eq, Show)

makeLenses ''ScoredPrimer
