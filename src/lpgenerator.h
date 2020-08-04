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
		lpGenerator::onlyone = onlyone;
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
	void predictProteinWithKNNCoil(char* AlphaIF, char* BetaIF,char* CoilIF, char* ProteinIF, std::string, int, int rand,bool firstIteration, int dirnum,bool ss,bool weights,bool dirr);
};

#endif // LPGENERATOR_H__»
