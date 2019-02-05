#!/bin/bash

lenwords=23
rand=1
numthreads=1
k=5
directorynum=4
firstitr=1
tag=output
seq=$1
python split.py $seq
iseq=$seq
seq=${seq}.seq
./lpgenerator PredictProteinKnnCoil $rand alpha.words beta.words other.words ${seq} $lenwords ${tag}/${seq}.knn $k $firstitr ${directorynum}
python getSecondaryStructure.py ${tag}/${seq}.knn ${tag}/${seq}.knnss 
python newpredictiontechniquecoil.py ${tag}/${seq}.knnss ${tag}/${seq}.probs $tag $k 
python convertprobabilities.py ${tag}/${seq}.probs $iseq ${tag}/${iseq}.ss
