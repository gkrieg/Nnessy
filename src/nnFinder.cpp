/* ********************************************************
 * Author: Spencer Krieger
 * Date: June 21st, 2016
 * For a description, please see the header file.
 * ********************************************************/

#include <iostream>
#include <fstream>
#include <string.h>
#include "nnFinder.h"
extern "C" {
	#include "dispersionTree.h"
	#include "metricspaces.h"
#include "word.h"
}
const char* residue = "ARNDCQEGHILKMFPSTWYVBZX-";


using namespace std;



MetricSpace *MS_QueriedWords;
MetricSpace *MS_QueryWords;
ifstream queriedWords;
ifstream queryWords;
int typeQuery;    /* 1 means querying for targets, 2 means querying for imposters */
int queryingClass;
int classToQuery;
int iteration;
string suffix;
string writeTo;



void nnFinder::makeDispersionTree(FILE *queriedWords, string type, int num, bool firstIteration, int rand){
	string input = "PROTEIN_WORDS";
        char* inputtype = new char[input.size() + 1];
        copy(input.begin(), input.end(), inputtype);
        inputtype[input.size()] = '\0';
	MetricSpace *MSQueriedWords = CreateMetricSpace(inputtype, queriedWords);
    string alphadfstring = "dfs/alphadis" + to_string(rand) + ".h";
    string betadfstring = "dfs/betadis" + to_string(rand) + ".h";
    string coildfstring = "dfs/coildis" + to_string(rand) + ".h";
    char* adf = new char[alphadfstring.size() + 1];
    copy(alphadfstring.begin(), alphadfstring.end(), adf);
    adf[alphadfstring.size()] = '\0';
    char* bdf = new char[betadfstring.size() + 1];
    copy(betadfstring.begin(), betadfstring.end(), bdf);
    bdf[betadfstring.size()] = '\0';
    char* cdf = new char[coildfstring.size() + 1];
    copy(coildfstring.begin(), coildfstring.end(), cdf);
    cdf[coildfstring.size()] = '\0';
    //alphadf = adf;
    //betadf = bdf;
    if (!firstIteration){
        setdfs(adf,bdf);
        int cset = setcoildf(cdf);
    }
	if (type == "alpha"){
		if (!firstIteration)
		    MSQueriedWords->dist = wordDist2Alpha;
        else
            MSQueriedWords->dist = WordDist;
		if (num == 1)
			AlphaDT = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 2)
			AlphaDT2 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 3)
			AlphaDT3 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
	}
	else if (type == "beta") {
		if (!firstIteration)
			MSQueriedWords->dist = wordDist2Beta;
        else
            MSQueriedWords->dist = WordDist;
		if (num == 1)
			BetaDT = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 2)
			BetaDT2 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 3)
			BetaDT3 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
	}
    else if (type == "alphaweighted"){
		if (!firstIteration)
		    MSQueriedWords->dist = WeightedAlphaDistance;
        else
            MSQueriedWords->dist = WeightedWordDist2;
		if (num == 1)
			AlphaDT = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 2)
			AlphaDT2 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 3)
			AlphaDT3 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
	}
	else if (type == "betaweighted") {
		if (!firstIteration)
			MSQueriedWords->dist = WeightedBetaDistance;
        else
            MSQueriedWords->dist = WeightedWordDist2;
		if (num == 1)
			BetaDT = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 2)
			BetaDT2 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 3)
			BetaDT3 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
	}
    else if (type == "alphacombined") {
        if (!firstIteration)
            MSQueriedWords->dist = SecondaryDist;
        else
            MSQueriedWords->dist = SecondaryDist;
		if (num == 1)
			AlphaDT = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 2)
			AlphaDT2 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 3)
			AlphaDT3 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
    }
    else if (type == "betacombined") {
        if (!firstIteration)
            MSQueriedWords->dist = SecondaryDist;
        else
            MSQueriedWords->dist = SecondaryDist;
		if (num == 1)
			BetaDT = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 2)
			BetaDT2 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 3)
			BetaDT3 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
    }
    else if (type == "coilcombined") {
        if (!firstIteration)
            MSQueriedWords->dist = SecondaryDist;
        else
            MSQueriedWords->dist = SecondaryDist;
        CoilDT = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
    }
	else if (type == "coil") {
		if (!firstIteration)
			MSQueriedWords->dist = wordDist2Coil;
        else
            MSQueriedWords->dist = WordDist;
        CoilDT = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
	}
	else if (type == "coilweighted") {
		if (!firstIteration)
			MSQueriedWords->dist = WeightedCoilDistance;
        else
            MSQueriedWords->dist = WeightedWordDist2;
        CoilDT = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
	}

	delete[] inputtype;
}
void nnFinder::loadDispersionTree(FILE *queriedWords, string filenamestring, string type, int rand){
	string input = "PROTEIN_WORDS";
	char* inputtype = new char[input.size() + 1];
	copy(input.begin(), input.end(), inputtype);
	inputtype[input.size()] = '\0';
	char* filename = new char[filenamestring.size() + 1];
	copy(filenamestring.begin(), filenamestring.end(), filename);
	filename[filenamestring.size()] = '\0';
	MetricSpace *MSQueriedWords = CreateMetricSpace(inputtype, queriedWords);
    string alphadfstring = "dfs/alphadis" + to_string(rand) + ".h";
    string betadfstring = "dfs/betadis" + to_string(rand) + ".h";
    string coildfstring = "dfs/coildis" + to_string(rand) + ".h";
    char* adf = new char[alphadfstring.size() + 1];
    copy(alphadfstring.begin(), alphadfstring.end(), adf);
    adf[alphadfstring.size()] = '\0';
    char* bdf = new char[betadfstring.size() + 1];
    copy(betadfstring.begin(), betadfstring.end(), bdf);
    bdf[betadfstring.size()] = '\0';
    char* cdf = new char[coildfstring.size() + 1];
    copy(coildfstring.begin(), coildfstring.end(), cdf);
    cdf[coildfstring.size()] = '\0';
    //alphadf = adf;
    //betadf = bdf;
    setdfs(adf,bdf);
    int cset = setcoildf(cdf);
	if (type == "alpha"){
		MSQueriedWords->dist = wordDist2Alpha;
		AlphaDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "beta"){
		MSQueriedWords->dist = wordDist2Beta;
		BetaDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "coil" && cset){
		MSQueriedWords->dist = wordDist2Coil;
		CoilDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "defaulta"){
		MSQueriedWords->dist = WordDist;
		AlphaDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "defaultb"){
		MSQueriedWords->dist = WordDist;
		BetaDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "defaultc"){
		MSQueriedWords->dist = WordDist;
		CoilDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "defaultaw"){
		MSQueriedWords->dist = WeightedWordDist2;
		AlphaDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "defaultbw"){
		MSQueriedWords->dist = WeightedWordDist2;
		BetaDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "defaultcw"){
		MSQueriedWords->dist = WeightedWordDist2;
		CoilDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "ass"){
		MSQueriedWords->dist = SecondaryDist2;
		AlphaDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "bss"){
		MSQueriedWords->dist = SecondaryDist2;
		BetaDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "css"){
		MSQueriedWords->dist = SecondaryDist2;
		CoilDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "alphaweighted"){
		MSQueriedWords->dist = WeightedAlphaDistance;
		AlphaDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "betaweighted"){
		MSQueriedWords->dist = WeightedBetaDistance;
		BetaDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "coilweighted" && cset){
		MSQueriedWords->dist = WeightedCoilDistance;
		CoilDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	delete[] inputtype;
	delete[] filename;

}
void nnFinder::destroyDTree(char type){
    if (type == 'A'){
        destroyDispersionTree(AlphaDT);
        AlphaDT = NULL;
    }
    else if (type == 'B'){
        destroyDispersionTree(BetaDT);
        BetaDT = NULL;
    }
    else{
        destroyDispersionTree(CoilDT);
        CoilDT = NULL;
    }
}
void nnFinder::saveDispersionTree(char* filename, string type){

	if(type == "alpha"){
		writeDispersionTreeToFile(filename, AlphaDT);
	}
	else if (type == "beta"){
		writeDispersionTreeToFile(filename, BetaDT);
	}
	else if (type == "coil"){
		writeDispersionTreeToFile(filename, CoilDT);
	}
}
void nnFinder::queryAlpha(int numNeighbors, MetricSpace *MS, Pointer*** AllNeighbors, int treenum){
	const int nump = MS->numpoints;
	(*AllNeighbors) = new Pointer*[nump];
	for (int numrows = 0; numrows < nump; numrows++){
		(*AllNeighbors)[numrows] = new Pointer[numNeighbors];
	}
	Pointer *words = MS->points;
	unsigned int i;
	int num;
	for ( i = 0; i < MS->numpoints; i++) {
		if (treenum == 1){
			num = knnSearch(AlphaDT, words[i], numNeighbors, (*AllNeighbors)[i]);
            //char* neighb = (char*) (AllNeighbors[i][0]);
            //printf("%s %s",(char*) words[i],neighb);
            }
		else if (treenum == 2)
			num = knnSearch(AlphaDT2, words[i], numNeighbors, (*AllNeighbors)[i]);
		else if (treenum == 3)
			num = knnSearch(AlphaDT3, words[i], numNeighbors, (*AllNeighbors)[i]);

	}

}
int nnFinder::queryAlphaRange(int sizeTree, MetricSpace *MS, Pointer*** AllNeighbors, int treenum){
	const int nump = MS->numpoints;
	(*AllNeighbors) = new Pointer*[nump];
	for (int numrows = 0; numrows < nump; numrows++){
		(*AllNeighbors)[numrows] = new Pointer[sizeTree];
	}
	Pointer *words = MS->points;
	unsigned int i;
	int num;
    char* nullstring = "null";
	for ( i = 0; i < MS->numpoints; i++) {
		if (treenum == 1){
			num = rangeSearch(AlphaDT, words[i], 5, (*AllNeighbors)[i]);
            (*AllNeighbors)[i][num] = nullstring;
            }
		else if (treenum == 2){
			num = rangeSearch(AlphaDT2, words[i], 5, (*AllNeighbors)[i]);
            (*AllNeighbors)[i][num] = nullstring;
    }
		else if (treenum == 3){
			num = rangeSearch(AlphaDT3, words[i], 5, (*AllNeighbors)[i]);
            (*AllNeighbors)[i][num] = nullstring;
    }

	}
    return num;

}
int nnFinder::queryBetaRange(int sizeTree, MetricSpace *MS, Pointer*** AllNeighbors, int treenum){
	const int nump = MS->numpoints;
	(*AllNeighbors) = new Pointer*[nump];
	for (int numrows = 0; numrows < nump; numrows++){
		(*AllNeighbors)[numrows] = new Pointer[sizeTree];
	}
    char* nullstring = "null";
	Pointer *words = MS->points;
	unsigned int i;
	int num;
	for ( i = 0; i < MS->numpoints; i++) {
		if (treenum == 1){
			num = rangeSearch(BetaDT, words[i], 5, (*AllNeighbors)[i]);
            (*AllNeighbors)[i][num] = nullstring;
            }
		else if (treenum == 2){
			num = rangeSearch(BetaDT2, words[i], 5, (*AllNeighbors)[i]);
            (*AllNeighbors)[i][num] = nullstring;
        }
		else if (treenum == 3){
			num = rangeSearch(BetaDT3, words[i], 5, (*AllNeighbors)[i]);
            (*AllNeighbors)[i][num] = nullstring;
        }

	}
    return num;

}
int nnFinder::queryCoilRange(int sizeTree, MetricSpace *MS, Pointer*** AllNeighbors, int treenum){
	const int nump = MS->numpoints;
    char* nullstring = "null";
	(*AllNeighbors) = new Pointer*[nump];
	for (int numrows = 0; numrows < nump; numrows++){
		(*AllNeighbors)[numrows] = new Pointer[sizeTree];
	}
	Pointer *words = MS->points;
	unsigned int i;
	int num;
	for ( i = 0; i < MS->numpoints; i++) {
		if (treenum == 1){
			num = rangeSearch(CoilDT, words[i], 5, (*AllNeighbors)[i]);
            (*AllNeighbors)[i][num] = nullstring;
            //char* neighb = (char*) (AllNeighbors[i][0]);
            //printf("%s %s",(char*) words[i],neighb);
            }
		else if (treenum == 2){
			num = rangeSearch(CoilDT, words[i], 5, (*AllNeighbors)[i]);
            (*AllNeighbors)[i][num] = nullstring;
        }
		else if (treenum == 3){
			num = rangeSearch(CoilDT, words[i], 5, (*AllNeighbors)[i]);
            (*AllNeighbors)[i][num] = nullstring;
        }

	}
    return num;

}
void nnFinder::queryCoil(int numNeighbors, MetricSpace *MS, Pointer*** AllNeighbors, int treenum){
	const int nump = MS->numpoints;
	(*AllNeighbors) = new Pointer*[nump];
	for (int numrows = 0; numrows < nump; numrows++){
		(*AllNeighbors)[numrows] = new Pointer[numNeighbors];
	}
	Pointer *words = MS->points;
	unsigned int i;
	int num;
	for ( i = 0; i < MS->numpoints; i++) {
        num = knnSearch(CoilDT, words[i], numNeighbors, (*AllNeighbors)[i]);

	}

}
void nnFinder::queryBeta(int numNeighbors, MetricSpace *MS, Pointer*** AllNeighbors, int treenum){
	const int nump = MS->numpoints;
	(*AllNeighbors) = new Pointer*[nump];
	for (int numrows = 0; numrows < nump; numrows++){
		(*AllNeighbors)[numrows] = new Pointer[numNeighbors];
	}
	Pointer *words = MS->points;
	unsigned int i;
	int num;
	for ( i = 0; i < MS->numpoints; i++) {
		if (treenum == 1)
			num = knnSearch(BetaDT, words[i], numNeighbors, (*AllNeighbors)[i]);
		else if (treenum == 2)
			num = knnSearch(BetaDT2, words[i], numNeighbors, (*AllNeighbors)[i]);
		else if (treenum == 3)
			num = knnSearch(BetaDT3, words[i], numNeighbors, (*AllNeighbors)[i]);

	}
}

