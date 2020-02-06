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
    case runExcept $ designOligsAA seed conf dna of
        Left err    -> trace ("error:" ++ err) $ return ()
        Right value -> trace ("value:" ++ prettyOligSet value) $ return ()
    return ()
