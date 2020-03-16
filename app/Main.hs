{-# LANGUAGE OverloadedStrings #-}

module Main where
import Bio.Tools.Sequence.OligoDesigner.Algo (designOligsAA, getRandomSeed)
import Control.Monad.Except (runExcept)
import Bio.Tools.Sequence.OligoDesigner.Types (OligsDesignerConfig(..))
import Bio.Tools.Sequence.OligoDesigner.Prettifier (prettyOligSet)
import Data.Default (def)
import Debug.Trace (trace)

--for performance testing
main :: IO ()
main = do
    seed <- getRandomSeed
    let conf = OligsDesignerConfig def def 0.7 0.3 0 1 1
    let dna = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDRTHTCPPCPAPELLGGPSVFLFPPKPKDTLYITREPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
    let result = runExcept $ designOligsAA seed conf dna
    case result of
        Left err    -> print ("error:" ++ err)
        Right value -> print ("value:" ++ prettyOligSet value)
