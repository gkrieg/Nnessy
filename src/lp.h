/**@file	/gsfs1/xdisk/skrieger/code/mgenerator/lpgenerator.h
 * @author	skrieger
 * @version	702
 * @date
 * 	Created:	Thu 23 Jun 2016 02:17:04 PM MST \n
 * 	Last Update:	Thu 23 Jun 2016 02:17:04 PM MST
 */
#ifndef LPGENERATOR_H_
#define LPGENERATOR_H_

#include <fstream>
#include "nnFinder.h"
#include <iomanip>
#include <iostream>
#include <string>
using namespace std;

const int ALPHABET_SIZE = 21;
const int numClasses = 3;
const double EPSILON = 0.0;
struct Mapping {
	int**** variableToIndex;
	int* indexToClass;
	int* indexToPosition;
	int* indexToFirstLetter;
	int* indexToSecondLetter;
};




class lpGenerator {
private:


public:
	Mapping* mapping;
	Pointer** AANeighbors;
	Pointer** ABNeighbors;
	Pointer** BBNeighbors;
	Pointer** BANeighbors;
	Pointer** ACNeighbors;
	Pointer** BCNeighbors;
	nnFinder nnfinder6;
	nnFinder nnfinder1;
	nnFinder nnfinder2;
	nnFinder nnfinder3;
	nnFinder nnfinder4;
	nnFinder nnfinder5;
	int k;
	int l;
	int TRSSize;
	int lengthOfWord;
	int numClasses;
	int numWordsAClass;
	int numWordsBClass;
	int numWordsCClass;
	bool onlyone;
	lpGenerator(int numClasses, int k, int l, int TRSSize, int wordlength, bool onlyone ){
		lpGenerator::numClasses = numClasses;
		lpGenerator::k = k;
		lpGenerator::l = l;
		lpGenerator::TRSSize = TRSSize;
		lpGenerator::lengthOfWord = wordlength;
        printf("wordlength is: %d",wordlength);
		lpGenerator::onlyone = onlyone;
		cout << numClasses << " " << k << " " << l << " " << TRSSize << endl;
		AANeighbors = NULL;
		ABNeighbors = NULL;
		BBNeighbors = NULL;
		ACNeighbors = NULL;
		BANeighbors = NULL;
		BCNeighbors = NULL;
	}

	~lpGenerator(){
		if (AANeighbors != NULL) {
			for (int i = 0;i < TRSSize;i++){
				delete[] AANeighbors[i];
			}
			delete[] AANeighbors;
		}
		if (!onlyone) {
			for (int i = 0;i < TRSSize;i++){
				delete[] ABNeighbors[i];
			}
			delete[] ABNeighbors;
		}
		if (!onlyone) {
			for (int i = 0;i < TRSSize;i++){
				delete[] BANeighbors[i];
			}
			delete[] BANeighbors;
		}
		if (!onlyone) {
			for (int i = 0;i < TRSSize;i++){
				delete[] BBNeighbors[i];
			}
			delete[] BBNeighbors;
		}
		if (!onlyone) {
			for (int i = 0;i < TRSSize;i++){
				delete[] ACNeighbors[i];
			}
			delete[] ACNeighbors;
		}
		if (!onlyone) {
			for (int i = 0;i < TRSSize;i++){
				delete[] BCNeighbors[i];
			}
			delete[] BCNeighbors;
		}
	}
	void findNeighbors(char* AlphaInputFile, char* BetaInputFile, char* AlphaTS, char* BetaTS, char* CoilTS, bool firstIteration, int rand);
    void printDFs();
	void createMapping();
	void writeObjectiveFunction(string alphaTSname, string betaTSname, string coilTSname, ofstream* output, int t_ssize);
	void writeSigmaFunctions(ofstream* output);
	void writeErrorFunctions(string alphaTSname, string betaTSname, string coilTSname, ofstream* output);
	void writeLP( string alphaTSname, string betaTSname, string coilTSname, ofstream* output);
	void deleteMapping();
    void findDistances(string QueryIF,string outputfile);
    void findDistances2(char*,char*, string ,int rand);
	void findNeighborsOneClass(char* TreeIF, char* QueryIF, int numneighbors, string outfile, bool firstIteration, int type, int rand);
	void findNeighborsOneClassCombined(char* TreeIF, char* QueryIF, int numneighbors, string outfile, bool firstIteration, int type, int rand, float epsilon, float rho);
	void findNeighborsCoverset(char* TreeIF, char* QueryIF, int numneighbors, string outfile, bool firstIteration, int type, int rand);
	void findNeighborsOneClassDistance(char* TreeIF, char* ,char* QueryIF, int numneighbors, string outfile, bool firstIteration, int type, int rand);
	void findNeighborsOneClassTwoDistance(char* AlphaIF, char* BetaIF, char* QueryIF, int numneighbors, string outfile, bool firstIteration, int type, int rand);
	void buildSaveTrees(char* alpha, char* beta, char* alphaSave, char* betaSave);
	void predict(char*, char*, char*, char, double, double, string, int rand);
	void predictProtein(char* AlphaIF, char* BetaIF, char* ProteinIF,  double, double, std::string, int, int rand);
	void predictProteinWithKNN(char* AlphaIF, char* BetaIF, char* ProteinIF,  double, double, std::string, int, int rand);
	void combineFiles();
	void writeErrorFunctionsFromFiles(ofstream* output);
	void threadalphafunc(int k, MetricSpace* MS, Pointer*** Neighbors, int treenum) {
		if (treenum == 1)
			nnfinder1.queryAlpha(k, MS, Neighbors, treenum);
		else if (treenum == 2)
			nnfinder2.queryAlpha(k, MS, Neighbors, treenum);
		else if (treenum == 3)
			nnfinder3.queryAlpha(k, MS, Neighbors, treenum);
		cout << "donethreadalphafunc" << endl;
	}
	void threadbetafunc(int k, MetricSpace* MS, Pointer*** Neighbors, int treenum) {
		if (treenum == 1)
			nnfinder4.queryBeta(k, MS, Neighbors, treenum);
		else if (treenum == 2)
			nnfinder5.queryBeta(k, MS, Neighbors, treenum);
		else if (treenum == 3)
			nnfinder6.queryBeta(k, MS, Neighbors, treenum);
		cout << "donethreadbetafunc" << endl;
	}
	void buildalphatree(FILE* alphaFile, int num, bool firstIteration,int rand) {
		if (num == 1)
			nnfinder1.makeDispersionTree(alphaFile, "alpha", num, firstIteration,rand);
		else if (num == 2)
			nnfinder2.makeDispersionTree(alphaFile, "alpha", num, firstIteration,rand);
		else if (num == 3)
			nnfinder3.makeDispersionTree(alphaFile, "alpha", num, firstIteration,rand);
		cout << "donebuildaketatree" << endl;
	}
	void buildcombinedalphatree(FILE* alphaFile, int num, bool firstIteration,int rand) {
		if (num == 1)
			nnfinder1.makeDispersionTree(alphaFile, "alphacombined", num, firstIteration, rand);
		else if (num == 2)
			nnfinder2.makeDispersionTree(alphaFile, "alphacombined", num, firstIteration, rand);
		else if (num == 3)
			nnfinder3.makeDispersionTree(alphaFile, "alphacombined", num, firstIteration, rand);
		cout << "donebuildaketatree" << endl;
	}
	void buildbetatree(FILE* betaFile, int num, bool firstIteration,int rand) {

		if (num == 1)
			nnfinder4.makeDispersionTree(betaFile, "beta", num, firstIteration, rand);
		else if (num == 2)
			nnfinder5.makeDispersionTree(betaFile, "beta", num, firstIteration, rand);
		else if (num == 3)
			nnfinder6.makeDispersionTree(betaFile, "beta", num, firstIteration, rand);
		cout << "donebuildbetatree" << endl;
	}
	void buildcombinedbetatree(FILE* betaFile, int num, bool firstIteration,int rand) {

		if (num == 1)
			nnfinder4.makeDispersionTree(betaFile, "betacombined", num, firstIteration, rand);
		else if (num == 2)
			nnfinder5.makeDispersionTree(betaFile, "betacombined", num, firstIteration, rand);
		else if (num == 3)
			nnfinder6.makeDispersionTree(betaFile, "betacombined", num, firstIteration, rand);
		cout << "donebuildbetatree" << endl;
	}
};

#endif // LPGENERATOR_H__»
