module Main where

import qualified Control.Applicative as CA
import qualified Data.Array.Repa as R
import Data.Attoparsec.Text as A
import GHC.Stack
import GPLVM.PPCA
import Prelude (read)
import System.IO hiding (putStrLn, readFile)
import System.Random
import Universum

main :: (HasCallStack) => IO ()
main = do
    [file, num] <- getArgs
    let num' = read num
    txt <- readFile file
    let (Right tParsedData) = A.parseOnly parseTestData txt
        parsedTestData = concat tParsedData
    gen <- getStdGen
    let matr = trace (show @String $ length parsedTestData) $ R.delay $ R.fromListUnboxed (R.Z R.:. 4 R.:. 150) parsedTestData
        ppca = makePPCA matr True 4(Left num') gen
        w = _W ppca
        lkh = _finalExpLikelihood ppca
    putStrLn @String $ show (R.computeS w :: R.Array R.U R.DIM2 Double)
    putStrLn @String $ show lkh


parseTestData :: Parser [[Double]]
parseTestData = transpose <$> sepBy (sepBy (A.double <|> (((read . toString) <$> A.string "NaN") :: A.Parser Double)) (CA.many $ char ' ')) (CA.many (char ' ') <* (string "\n" CA.<|> string "\r\n"))
