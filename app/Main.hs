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
    let matr = R.delay $ R.fromListUnboxed (R.Z R.:. 4 R.:. 150) parsedTestData
        pca = makePPCA matr 4 Nothing (Just num') gen
        likelihood = _finalExpLikelihood pca
    putStrLn @String $ show likelihood

parseTestData :: Parser [[Double]]
parseTestData = transpose <$> sepBy (sepBy A.double (CA.many $ char ' ')) (CA.many (char ' ') <* (string "\n" CA.<|> string "\r\n"))
