#!/bin/bash

touchstonef=$1
lengthofword=23          # length of words to be generated
python3 filegen.py $touchstonef $lengthofword
#python3 filegen8.py $touchstonef $lengthofword
python removeduplicates.py
python countNumberOfLinesInFile.py alpha.words beta.words other.words
