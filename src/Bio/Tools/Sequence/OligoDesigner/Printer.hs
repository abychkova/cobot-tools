module Bio.Tools.Sequence.OligoDesigner.Printer (
    printResult
    ,buildStr
    ,buildStr'
) where

import Bio.Tools.Sequence.OligoDesigner.Types (OligsDesignerConfig(..), OligSet(..), Olig(..))
import Text.Regex.TDFA (Regex, makeRegex)
import Bio.Tools.Sequence.CodonOptimization (gcContentDesired)
import Bio.Tools.Sequence.OligoDesigner.Scorer (commonScore, gcContent, gcContentScoreByOligs, rnaScore, oligsGCContentDifference)
import Bio.Tools.Sequence.OligoDesigner.Utils (assemble)
import Bio.NucleicAcid.Nucleotide (DNA(..))
import Bio.Tools.Sequence.CodonOptimization.Algo (scoreByWindow)
import Bio.Tools.Sequence.OligoDesigner.Prettifier (prettyDNA)
import Bio.Tools.Sequence.CodonOptimization.Types (defaultForbiddenRegexp)
import Debug.Trace (trace)


printResult :: OligsDesignerConfig -> String -> OligSet -> IO ()  
printResult conf app ols = print $ buildStr conf app ols

buildStr :: OligsDesignerConfig -> String -> OligSet -> String
buildStr (OligsDesignerConfig codonConf _ rnaF oligsGCF gcF _ _) app ols= str
  where
    target = gcContentDesired codonConf
    str = "," ++ app ++ "," ++ 
        show (ylabScore gcScoreVal rnaScoreVal) ++ "," ++
        show (commonScore' gcScoreVal rnaScoreVal oligsGCValue) ++ "," ++
        show rnaScoreVal ++ "," ++ 
        show gcScoreVal ++ "," ++ 
        show oligsGCValue ++ "," ++ 
        show (gcContent ols) ++ "," ++ 
        show (tempDifference ols) ++ "," ++ 
        show (scoreByWindow codonConf regexes (assemble ols)) ++ "," ++ 
        prettyDNA(assemble ols)
            
    commonScore' :: Double -> Double -> Double -> Double
    commonScore' gcScore rnaScore oligsGC = rnaScore * gcScore / oligsGC
      where
        oligsGCValue = oligsGCContentDifference ols    
           
    regexes :: [Regex]
    regexes = map makeRegex defaultForbiddenRegexp :: [Regex]
       
    gcScore :: Double -> Double        
    gcScore = gcContentScoreByOligs ols
    
    rnaScoreVal = realToFrac $ rnaScore ols
    gcScoreVal = gcScore target
    oligsGCValue = oligsGCContentDifference ols    
            
    ylabScore :: Double -> Double -> Double
    ylabScore gcScore rnaScore = if rnaScore > 0 then gcScore ** 50 * rnaScore else rnaScore
            
    tempDifference :: OligSet -> Double
    tempDifference (OligSet fwd rvsd _) = maximum temps - minimum temps
      where
        oligs = fwd ++ rvsd
        temps = map temp oligs
            
    temp :: Olig -> Double
    temp (Olig dna _ _) = 64.9 + 41 * (gc - 16.4) / len
      where
        len = realToFrac $ length dna
        gc = realToFrac $ length $ filter (\nk -> nk == DC || nk == DG) dna
    
    gcContent :: OligSet -> Double
    gcContent oligs = gc / realToFrac (length dna)
      where
        dna = assemble oligs
        gc = realToFrac $ length $ filter (\nk -> nk == DC || nk == DG) dna
        
buildStr' :: OligsDesignerConfig -> String -> OligSet -> String
buildStr' conf@(OligsDesignerConfig codonConf _ rnaF oligsGCF gcF _ _) app ols= str
  where
    target = gcContentDesired codonConf
    str = "," ++ app ++ "," ++ 
        show (ylabScore gcScoreVal rnaScoreVal) ++ "," ++
        show (commonScore' gcScoreVal rnaScoreVal oligsGCValue) ++ "," ++
--        show (commonScore conf ols) ++ "," ++
        show rnaScoreVal ++ "," ++ 
        show gcScoreVal ++ "," ++ 
        show oligsGCValue ++ "," ++ 
        show (gcContent ols) ++ "," ++ 
        show (tempDifference ols) ++ "," ++ 
        show (scoreByWindow codonConf regexes (assemble ols))
            
    commonScore' :: Double -> Double -> Double -> Double
    commonScore' gcScore rnaScore oligsGC = rnaScore * gcScore / oligsGC 
           
    regexes :: [Regex]
    regexes = map makeRegex defaultForbiddenRegexp :: [Regex]
       
    gcScore :: Double -> Double        
    gcScore = gcContentScoreByOligs ols
    
    rnaScoreVal = realToFrac $ rnaScore ols
    gcScoreVal = gcScore target
    oligsGCValue = oligsGCContentDifference ols    
            
    ylabScore :: Double -> Double -> Double
    ylabScore gcScore rnaScore = if rnaScore > 0 then gcScore ** 50 * rnaScore else rnaScore
            
    tempDifference :: OligSet -> Double
    tempDifference (OligSet fwd rvsd _) = maximum temps - minimum temps
      where
        oligs = fwd ++ rvsd
        temps = map temp oligs
            
    temp :: Olig -> Double
    temp (Olig dna _ _) = 64.9 + 41 * (gc - 16.4) / len
      where
        len = realToFrac $ length dna
        gc = realToFrac $ length $ filter (\nk -> nk == DC || nk == DG) dna
    
    gcContent :: OligSet -> Double
    gcContent oligs = gc / realToFrac (length dna)
      where
        dna = assemble oligs
        gc = realToFrac $ length $ filter (\nk -> nk == DC || nk == DG) dna