module Bio.Tools.Sequence.Annotation.Algorithm where

import           Bio.NucleicAcid.Nucleotide           (Complementary (..),
                                                       DNA (..), symbol)
import Data.Text (Text)
import Bio.Chain.Alignment
import Bio.Chain.Alignment.Scoring
import Control.Arrow
import Data.Ord
import Data.List
import Debug.Trace (trace)
import           System.IO                        (BufferMode (..), IOMode (..),
                                                   hPutStrLn, hSetBuffering,
                                                   stdout, withFile)


mkHappy :: String -> [NamedA]
mkHappy s = map (second (align parameters s)) abSequences
  where
    parameters = SemiglobalAlignment blosum62 (AffineGap (-11) (-1))
    

type Named a = (Text, a)
type NamedA = Named (AlignmentResult String String)

abSequences :: [Named String]
abSequences = 
    [ ("IGHG1_HUMAN_CH1", "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV")
    , ("IGHG2_HUMAN_CH1", "ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSNFGTQTYTCNVDHKPSNTKVDKTV")
    , ("IGHG3_HUMAN_CH1", "ASTKGPSVFPLAPCSRSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYTCNVNHKPSNTKVDKRV")
    , ("IGHG4_HUMAN_CH1", "ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTKTYTCNVDHKPSNTKVDKRV")
    , ("IGHG1_HUMAN_Hinge", "EPKSCDKTHTCP")
    , ("IGHG2_HUMAN_Hinge", "ERKCCVECPPCP")
    , ("IGHG3_HUMAN_Hinge", "ELKTPLGDTTHTCPRCPEPKSCDTPPPCPRCPEPKSCDTPPPCPRCPEPKSCDTPPPCPRCP")
    , ("IGHG4_HUMAN_Hinge", "ESKYGPPCPSCP")
    , ("IGHG1_HUMAN_CH2", "PCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAK")
    , ("IGHG2_HUMAN_CH2", "APPVAGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTFRVVSVLTVVHQDWLNGKEYKCKVSNKGLPAPIEKTISKTK")
    , ("IGHG3_HUMAN_CH2", "APELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVQFKWYVDGVEVHNAKTKPREEQYNSTFRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKTK")
    , ("IGHG4_HUMAN_CH2", "APEFLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKGLPSSIEKTISKAK")
    , ("IGHG1_HUMAN_CH3", "GQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK")
    , ("IGHG2_HUMAN_CH3", "GQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDISVEWESNGQPENNYKTTPPMLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK")
    , ("IGHG3_HUMAN_CH3", "GQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESSGQPENNYNTTPPMLDSDGSFFLYSKLTVDKSRWQQGNIFSCSVMHEALHNRFTQKSLSLSPG")
    , ("IGHG4_HUMAN_CH3", "GQPREPQVYTLPPSQEEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSRLTVDKSRWQEGNVFSCSVMHEALHNHYTQKSLSLSLGK")
    , ("IGLC1", "GQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADGSPVKAGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS")
    , ("IGLC2", "GQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS")
    , ("IGLC2", "GQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHKSYSCQVTHEGSTVEKTVAPTECS")
    , ("IGKC", "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC")
    ]

{-
findBreakPoint :: String -> ((String, String), Float)
findBreakPoint s = (flip splitAt s . length . takeWhile isDelete . alignment $ withMaxScore, myScore withMaxScore)
  where
    result = snd <$> mkHappy s
    scores = (\aRes -> (fromIntegral . score $ aRes) / (fromIntegral . length . sequence2 $ aRes)) <$> result
    withMaxScore = maximumBy (comparing myScore) result 
-}

findBestAlignment :: String -> NamedA -- AlignmentResult String String --NamedA
findBestAlignment s = withMaxScore
  where
    result = mkHappy s
    withMaxScore = maximumBy (comparing (myScore . snd)) result 

mkMeHappy :: String -> [Named String]
mkMeHappy "" = []
mkMeHappy s = do
    let best = findBestAlignment s
    if myScore (snd best) > 1
      then do
        let (a, b, c) = splitSequence best
        concat [mkMeHappy a, [b], mkMeHappy c]
      else []

myScore :: AlignmentResult String String -> Float
myScore ar = (fromIntegral . score $ ar) / (fromIntegral . length . sequence2 $ ar)

isDelete :: Operation i j -> Bool
isDelete DELETE {} = True
isDelete _         = False

findStart :: AlignmentResult String String -> Int
findStart = length . takeWhile isDelete . alignment

findEnd :: AlignmentResult String String -> Int
findEnd ar = (length . sequence1 $ ar) - (length . takeWhile isDelete . reverse . alignment $ ar)

splitSequence :: NamedA -> (String, Named String, String)
splitSequence (name, ar) = (before, (name, current), after)
  where
    seq'   = sequence1 ar

    start' = findStart ar
    end'   = findEnd ar

    before  = take start' seq'
    current = take (end' - start') . drop start' $ seq'
    after   = drop end' seq'

---
heavy = "EVQLVQSGAEVKKPGASVKVSCKASGYTFTNYYIHWVRQAPGQGLEWIGWIYPGDGNTKYNEKEKGRATLTADTSTSTAYLELSSLRSEDTAVYYCARDSYSNYYFDYWGQGTLVTVSSASTKGPSVEPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVELFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKENWYVDGVEVHNAKTKPREEQYGSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLSCAVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLVSKLTVDKSRWQQGNVESCSVMHEALHNHYTQKSLSLSPGK" :: String

hinge = "EPKSCDKTHTCP" :: String

res = mkHappy heavy
x = snd $ res !! 4
