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
#include <map>
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
	Pointer** CANeighbors;
	Pointer** CCNeighbors;
	Pointer** CBNeighbors;
	nnFinder nnfinder6;
	nnFinder nnfinder7;
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
    map<string,int> countmap;
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
		CCNeighbors = NULL;
		CBNeighbors = NULL;
		CANeighbors = NULL;
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
		if (!onlyone) {
			for (int i = 0;i < TRSSize;i++){
				delete[] CCNeighbors[i];
			}
			delete[] CCNeighbors;
		}
		if (!onlyone) {
			for (int i = 0;i < TRSSize;i++){
				delete[] CANeighbors[i];
			}
			delete[] CANeighbors;
		}
		if (!onlyone) {
			for (int i = 0;i < TRSSize;i++){
				delete[] CBNeighbors[i];
			}
			delete[] CBNeighbors;
		}
	}
	void findNeighbors(char* AlphaInputFile, char* BetaInputFile, char* AlphaTS, char* BetaTS, char* CoilTS, bool firstIteration, int rand);
    void testDF(int,char*,char*);
    void printDFs();
	void createMapping();
	void writeObjectiveFunction(string alphaTSname, string betaTSname, string coilTSname, ofstream* output, int t_ssize);
	void writeSigmaFunctions(ofstream* output);
	void writeErrorFunctions(string alphaTSname, string betaTSname, string coilTSname, ofstream* output);
	void writeLP( string alphaTSname, string betaTSname, string coilTSname, ofstream* output);
    void buildmap(string);
	void deleteMapping();
    void findDistances(string QueryIF,string outputfile);
    void findDistances2(char*,char*, string ,int rand);
    void findCombinedDistances2(char*,char*, string ,int rand, float, float, float);
	void findNeighborsOneClass(char* TreeIF, char* QueryIF, int numneighbors, string outfile, bool firstIteration, int type, int rand, float rho, bool counts, bool load);
	void findNeighborsOneClassCombined(char* TreeIF, char* QueryIF, int numneighbors, string outfile, bool firstIteration, int type, int rand, float epsilon, float aarho, float ssrho);
	void findNeighborsCoverset(char* TreeIF, char* QueryIF, int numneighbors, string outfile, bool firstIteration, int type, int rand);
	void findNeighborsOneClassDistance(char* TreeIF, char* ,char* QueryIF, int numneighbors, string outfile, bool firstIteration, int type, int rand);
	void findNeighborsOneClassTwoDistance(char* AlphaIF, char* BetaIF, char* QueryIF, int numneighbors, string outfile, bool firstIteration, int type, int rand, bool load, bool weights, bool training);
	void findNeighborsCoilClassTwoDistance(char* AlphaIF, char* QueryIF, int numneighbors, string outfile, bool firstIteration, int type, int rand, bool load, bool weights, bool training);
	void findNeighborsOneClassTwoDistanceWeighted(char* AlphaIF, char* BetaIF, char* QueryIF, int numneighbors, string outfile, bool firstIteration, int type, int rand, float rho);
void findNeighborsOneClassTwoDistanceAve(char* AlphaIF, char* BetaIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand, bool counts);
void findNeighborsOneClassTwoDistanceAveLoad(char* AlphaIF, char* BetaIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand, bool counts);
void findNeighborsOneClassTwoDistanceAveLoadCoil(char* AlphaIF, char* BetaIF,char* CoilIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand, bool counts);
    void setshift(float s);
	void buildSaveTrees(char* alpha, char* beta, char* alphaSave, char* betaSave, int rand,bool weights);
	void buildSaveTreesCoil(char* alpha, char* beta,char* coil, char* alphaSave, char* betaSave,char* coilSave, int rand,bool weights);
	void buildSStrees(char* alpha, char* beta,char* coil, char* alphaSave, char* betaSave,char* coilSave, int rand,bool weights);
	void buildSaveTreesDefault(char* alpha, char* beta, char* alphaSave, char* betaSave, int rand,bool weights);
	void buildSaveTreesCoilDefault(char* alpha, char* beta,char* coil, char* alphaSave, char* betaSave,char* coilSave, int rand,bool weights);
	void predict(char*, char*, char*, char, double, double, string, int rand);
	void predictProtein(char* AlphaIF, char* BetaIF, char* ProteinIF,  double, double, std::string, int, int rand, bool firstIteration);
	void predictProteinLoadTree(char* AlphaIF, char* BetaIF, char* ProteinIF,  double, double, std::string, int, int rand, std::string Treefile, std::string Treefile2, bool parallel, bool weights);
	void predictProteinLoadTreeDefault(char* AlphaIF, char* BetaIF, char* ProteinIF,  double, double, std::string, int, int rand, std::string Treefile1,std::string Treefile2, bool parallel, bool weights, bool training);
	void predictProteinLoadTreeDefaultCoil(char* AlphaIF, char* BetaIF,char* CoilIF, char* ProteinIF,  double, double, std::string, int, int rand, std::string Treefile1,std::string Treefile2,std::string Treefile3, bool parallel, bool weights, bool training);
	void predictProteinLoadTreeCoil(char* AlphaIF, char* BetaIF,char* CoilIF, char* ProteinIF,  double, double, std::string, int, int rand, std::string Treefile1,std::string Treefile2,std::string Treefile3, bool parallel, bool weights, bool training);
	void predictProteinWithKNN(char* AlphaIF, char* BetaIF, char* ProteinIF,  double, double, std::string, int, int rand);
	void predictProteinWithKNNCoil(char* AlphaIF, char* BetaIF,char* CoilIF, char* ProteinIF, std::string, int, int rand,bool firstIteration, int dirnum,bool ss,bool weights,bool dirr);
    void gettargetsdist(int);
	void combineFiles();
	void writeErrorFunctionsFromFiles(ofstream* output);
    void getnewtargets(int);
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
	/*void* threadalphafunc2(void* tdp) {
        threaddata td = (struct threaddata) tdp;
		if (td.treenum == 1)
			nnfinder1.queryAlpha(td.numneighbors, td.MSAlpha, td.writespace, td.treenum);
		else if (td.treenum == 2)
			nnfinder2.queryAlpha(td.numneighbors, td.MSAlpha, td.writespace, td.treenum);
		else if (td.treenum == 3)
			nnfinder3.queryAlpha(td.numneighbors, td.MSAlpha, td.writespace, td.treenum);
		cout << "donethreadalphafunc" << endl;
	}
    
	void threadbetafunc2(int k, MetricSpace* MS, Pointer*** Neighbors, int treenum) {
		if (treenum == 1)
			nnfinder4.queryBeta(k, MS, Neighbors, treenum);
		else if (treenum == 2)
			nnfinder5.queryBeta(k, MS, Neighbors, treenum);
		else if (treenum == 3)
			nnfinder6.queryBeta(k, MS, Neighbors, treenum);
		cout << "donethreadbetafunc" << endl;
	}
    */
	void buildalphatree(FILE* alphaFile, int num, bool firstIteration,int rand) {
		if (num == 1)
			nnfinder1.makeDispersionTree(alphaFile, "alpha", num, firstIteration,rand);
		else if (num == 2)
			nnfinder2.makeDispersionTree(alphaFile, "alpha", num, firstIteration,rand);
		else if (num == 3)
			nnfinder3.makeDispersionTree(alphaFile, "alpha", num, firstIteration,rand);
		cout << "donebuildaketatree" << endl;
	}
	void buildalphatreeweighted(FILE* alphaFile, int num, bool firstIteration,int rand) {
		if (num == 1)
			nnfinder1.makeDispersionTree(alphaFile, "alphaweighted", num, firstIteration,rand);
		else if (num == 2)
			nnfinder2.makeDispersionTree(alphaFile, "alphaweighted", num, firstIteration,rand);
		else if (num == 3)
			nnfinder3.makeDispersionTree(alphaFile, "alphaweighted", num, firstIteration,rand);
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
	void buildbetatreeweighted(FILE* betaFile, int num, bool firstIteration,int rand) {

		if (num == 1)
			nnfinder4.makeDispersionTree(betaFile, "betaweighted", num, firstIteration, rand);
		else if (num == 2)
			nnfinder5.makeDispersionTree(betaFile, "betaweighted", num, firstIteration, rand);
		else if (num == 3)
			nnfinder6.makeDispersionTree(betaFile, "betaweighted", num, firstIteration, rand);
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
	void buildcoiltree(FILE* coilFile, int num, bool firstIteration,int rand) {
        nnfinder7.makeDispersionTree(coilFile, "coil", num, firstIteration,rand);
		cout << "donebuildcoilatree" << endl;
	}
	void buildcoiltreeweighted(FILE* coilFile, int num, bool firstIteration,int rand) {
        nnfinder7.makeDispersionTree(coilFile, "coilweighted", num, firstIteration,rand);
		cout << "donebuilcoilwtree" << endl;
	}
	void buildcombinedcoiltree(FILE* coilFile, int num, bool firstIteration,int rand) {
        nnfinder7.makeDispersionTree(coilFile, "coilcombined", num, firstIteration,rand);
		cout << "donebuilcoilwtree" << endl;
	}
};

#endif // LPGENERATOR_H__»
