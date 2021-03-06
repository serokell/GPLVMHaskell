module Main where

import qualified Control.Applicative as CA
import qualified Data.Array.Repa as R
import Data.Attoparsec.Text as A
import qualified GPLVM.PPCA as U
import GPLVM.TypeSafe.PPCA
import GPLVM.TypeSafe.Types
import Prelude (read)
import System.Environment
import System.IO hiding (putStrLn, readFile)
import System.Random
import Universum

main :: IO ()
main = do
    [file, num', ts] <- getArgs
    let num = read num'
    let typeSafe = read ts
    txt <- readFile file
    let (Right tParsedData) = A.parseOnly parseTestData txt
        parsedTestData = concat tParsedData
    gen <- getStdGen
    let matr = R.delay $ R.fromListUnboxed (R.Z R.:. 4 R.:. 150) parsedTestData
        ppca = makePPCATypeSafe matr 3 (Left num) gen
        uppca = U.makePPCA matr 3 (Left num) gen
        w = _W ppca
        wU = U._W uppca
        lkh = _finalExpLikelihood ppca
        lkhU = U._finalExpLikelihood uppca
        restored = _restoredMatrix ppca
        restoredU = U._restoredMatrix uppca
    if typeSafe
    then do
        putStrLn @String $ show (R.computeS w :: R.Array R.U R.DIM2 Double)
        putStrLn @String $ show lkh
        putStrLn @String $ show @String (R.computeS <$> restored :: Maybe (R.Array R.U R.DIM2 Double))
    else do
        putStrLn @String $ show (R.computeS wU :: R.Array R.U R.DIM2 Double)
        putStrLn @String $ show lkhU
        putStrLn @String $ show @String (R.computeS <$> restoredU :: Maybe (R.Array R.U R.DIM2 Double))

parseTestData :: Parser [[Double]]
parseTestData = transpose <$> sepBy (sepBy (A.double <|> (((read . toString) <$> A.string "NaN") :: A.Parser Double)) (CA.many $ char ' ')) (CA.many (char ' ') <* (string "\n" CA.<|> string "\r\n"))
