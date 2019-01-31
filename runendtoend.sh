#!/bin/bash

lenwords=23
rand=1
PBS_ARRAY_INDEX=0
numthreads=1
#numthreads=4
k=5
directorynum=4
firstitr=1
tag=output
seq=$1
###Need to first do some splitting based on parameters
#python splitfile.py $seq $numthreads
###Then we run the neighbor finding stuff
#for i in `seq 0 $numthreads`; do
#./lpgenerator PredictProteinKnnCoil $rand alpha.words beta.words other.words ${seq}_${numthreads}_${i} $lenwords predictions/${rand}${tag}/predictknn${seqnum}.txt_${i} $k $firstitr ${directorynum}>> stuff.txt &
python split.py $seq
seq=${seq}.seq
./lpgenerator PredictProteinKnnCoil $rand alpha.words beta.words other.words ${seq} $lenwords ${tag}/${seq}.knn $k $firstitr ${directorynum}
#done
#wait
###Then we find the secondary structure for that specific sequence
#for i in `seq 0 $numthreads`; do
python getSecondaryStructure.py ${tag}/${seq}.knn ${tag}/${seq}.knnss $directorynum 
#done
#wait
###Then we do the newpred
#for i in `seq 0 $numthreads`; do
python newpredictiontechniquecoil.py ${tag}/${seq}.knnss ${tag}/${seq}.probs $tag $k 
#done
#wait
#python combinefiles.py predictions/${rand}${tag}/predictknn${seqnum}.txt $numthreads
#python combinefiles.py predictions/${rand}${tag}/predictknn${seqnum}.txtss $numthreads
#python combinefiles.py predictions/${rand}${tag}/probabilities${seqnum} $numthreads
#TARGET=$(date)
#one=$(diff)
#echo $seq $numthreads $one >> times2.txt
###numthreads
