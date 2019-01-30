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
//#include "alphasubstitutionScores.h"
//#include "betasubstitutionScores.h"
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
    cout << "printing out the word" << (char*) MSQueriedWords->points[2];
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
    cout << endl << "setting distance functions with rand = " << rand << endl;
    if (!firstIteration){
        setdfs(adf,bdf);
        int cset = setcoildf(cdf);
    }
    cout << "done setting distance functions" << endl;
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
            MSQueriedWords->dist = WeightedWordDist;
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
            MSQueriedWords->dist = WeightedWordDist;
		if (num == 1)
			BetaDT = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 2)
			BetaDT2 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 3)
			BetaDT3 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
	}
    else if (type == "alphacombined") {
        if (!firstIteration)
            MSQueriedWords->dist = CombinedAlphaDistance;
        else
            MSQueriedWords->dist = CombinedWordDist;
		if (num == 1)
			AlphaDT = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 2)
			AlphaDT2 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 3)
			AlphaDT3 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
    }
    else if (type == "betacombined") {
        if (!firstIteration)
            MSQueriedWords->dist = CombinedBetaDistance;
        else
            MSQueriedWords->dist = CombinedWordDist;
		if (num == 1)
			BetaDT = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 2)
			BetaDT2 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
		else if (num == 3)
			BetaDT3 = createDispersionTree(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, 2, 0);
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
            MSQueriedWords->dist = WeightedWordDist;
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
    cout << endl << "setting distance functions with rand = " << rand << endl;
    setdfs(adf,bdf);
    int cset = setcoildf(cdf);
    cout << endl << "cset " << cset << endl;
    cout << "done setting distance functions" << endl;
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
		MSQueriedWords->dist = WeightedWordDist;
		AlphaDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "defaultbw"){
		MSQueriedWords->dist = WeightedWordDist;
		BetaDT = readDispersionTreeFromFile(MSQueriedWords->points, MSQueriedWords->dist, MSQueriedWords->numpoints, filename);
	}
	else if (type == "defaultcw"){
		MSQueriedWords->dist = WeightedWordDist;
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
	cout << "numpoints alpha is " << nump << endl;
	(*AllNeighbors) = new Pointer*[nump];
	for (int numrows = 0; numrows < nump; numrows++){
		(*AllNeighbors)[numrows] = new Pointer[numNeighbors];
	}
	Pointer *words = MS->points;
	unsigned int i;
	int num;
	for ( i = 0; i < MS->numpoints; i++) {
		if (treenum == 1)
			num = knnSearch(AlphaDT, words[i], numNeighbors, (*AllNeighbors)[i]);
		else if (treenum == 2)
			num = knnSearch(AlphaDT2, words[i], numNeighbors, (*AllNeighbors)[i]);
		else if (treenum == 3)
			num = knnSearch(AlphaDT3, words[i], numNeighbors, (*AllNeighbors)[i]);

	}
	cout << i << endl << num << endl;

}
void nnFinder::queryCoil(int numNeighbors, MetricSpace *MS, Pointer*** AllNeighbors, int treenum){
	const int nump = MS->numpoints;
	cout << "numpoints coil is " << nump << endl;
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
	cout << i << endl << num << endl;

}
void nnFinder::queryBeta(int numNeighbors, MetricSpace *MS, Pointer*** AllNeighbors, int treenum){
	const int nump = MS->numpoints;
	cout << "numpoints is " << nump << endl;
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
	cout << i << endl << num << endl;
}

/*void nnFinder::printDFs(){
    int length = 22;
    cout << "a1{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a1[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a2{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a2[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a3{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a3[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a4{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a4[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a5{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a5[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a6{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a6[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a7{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a7[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a8{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a8[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a9{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a9[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a10{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a10[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a11{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a11[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a12{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a12[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a13{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a13[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a14{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a14[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a15{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a15[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a16{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a16[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a17{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a17[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a18{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a18[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a19{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a19[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a20{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a20[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a21{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a21[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a22{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a22[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "a23{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << a23[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl << endl;
    cout << "bset now:" << endl;
    cout << "b1{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b1[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b2{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b2[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b3{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b3[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b4{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b4[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b5{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b5[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b6{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b6[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b7{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b7[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b8{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b8[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b9{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b9[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b10{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b10[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b11{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b11[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b12{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b12[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b13{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b13[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b14{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b14[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b15{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b15[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b16{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b16[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b17{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b17[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b18{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b18[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b19{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b19[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b20{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b20[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b21{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b21[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b22{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b22[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    cout << "b23{";
    for (int i = 0; i < length; i++){
        cout << residue[i] << " ";
        for (int j = 0; j < length; j++){
            cout << residue[j] << " ";
            cout << b23[i][j] << ",";
        }
        cout << endl;
    }
    cout << "}" << endl;
    printepsilon();
}
*/
