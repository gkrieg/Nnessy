/**@file	/gsfs1/xdisk/skrieger/code/mgenerator/main.cpp
 * @author	skrieger
 * @version	702
 * @date
 * 	Created:	Mon 27 Jun 2016 12:56:18 PM MST \n
 * 	Last Update:	Mon 27 Jun 2016 12:56:18 PM MST
 */

/*===========================================================================*/
/*===============================[ main for lpgen ]===============================*/
/*===========================================================================*/
#include "lpgenerator.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
using namespace std;
/*===========================================================================*/
/*===============================[ main for lpgen ]===============================*/
/*===========================================================================*/

int main(int argc, char* argv[]) {
    int rand = 0;
	if (argc == 12) {
		ofstream output;
		string filename = "words.lp";
		output.open(filename.c_str());
		int numClasses = atoi( argv[6]);
		int k = atoi(argv[7]);
		int l = atoi(argv[8]);
		int TRSSize = atoi(argv[9]);
		int wordlength = atoi(argv[10]);
		int firstIteration = atoi(argv[11]);
		bool isFirstIteration;
		if (firstIteration == 1)
			isFirstIteration = true;
		else
			isFirstIteration = false;
		lpGenerator lpgen = lpGenerator(numClasses, k, l, TRSSize, wordlength, false);
		char* alpha = argv[1];
		char* beta = argv[2];
		char* alphats = argv[3];
		char* betats = argv[4];
		char* coilts = argv[5];
		cout << "starting to find neighbors" << endl;
		lpgen.findNeighbors(alpha, beta, alphats, betats, coilts, isFirstIteration, rand);
		cout << "found all the neighbors, now writing to file" << endl;
		lpgen.writeLP(alphats, betats, coilts, &output);
		cout << "I'm at the end of lpgen" << endl;
		output.close();
	}
    //From manyjobs.sh
	else if (argc == 9){
		char* Tree = argv[1];
		char* Query = argv[2];
		int numneighbors = atoi(argv[3]);
		string outfile = argv[4];
		int numclasses = atoi(argv[5]);
		int wordlength = atoi(argv[6]);
		int firstIteration = atoi(argv[7]);
		bool isFirstIteration;
        int type = atoi(argv[8]);
		if (firstIteration == 1)
			isFirstIteration = true;
		else
			isFirstIteration = false;
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClass(Tree, Query, numneighbors, outfile, isFirstIteration,type, rand);
	}
	else if (argc == 15){
        cout << "I'm here in the find combined" << endl;
		char* Tree = argv[1];
		char* Query = argv[2];
		int numneighbors = atoi(argv[3]);
		string outfile = argv[4];
		int numclasses = atoi(argv[5]);
		int wordlength = atoi(argv[6]);
		int firstIteration = atoi(argv[7]);
		bool isFirstIteration;
        int type = atoi(argv[8]);
		if (firstIteration == 1)
			isFirstIteration = true;
		else
			isFirstIteration = false;
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassCombined(Tree, Query, numneighbors, outfile, isFirstIteration,type, rand);
	}
	else if (argc == 13){
		char* Tree = argv[1];
		char* Query = argv[2];
		int numneighbors = atoi(argv[3]);
		string outfile = argv[4];
		int numclasses = atoi(argv[5]);
		int wordlength = atoi(argv[6]);
		int firstIteration = atoi(argv[7]);
		bool isFirstIteration;
        int type = atoi(argv[8]);
        if (firstIteration == 1)
            isFirstIteration = true;
        else
            isFirstIteration = false;
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsCoverset(Tree, Query, numneighbors, outfile, isFirstIteration,type, rand);
	}
    //From twodist.sh
	else if (argc == 7){
        cout << "argc was 7" << endl;
		char* Alpha = argv[1];
		char* Beta = argv[2];
		char* Query = argv[3];
		int numneighbors = 2;
		string outfile = argv[4];
		int type = atoi(argv[5]);
		int wordlength = atoi(argv[6]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassTwoDistance(Alpha, Beta, Query, numneighbors, outfile, false, type, rand);
	}
	else if (argc == 8){
        cout << "argc was 8" << endl;
		char* Alpha = argv[1];
		char* Beta = argv[2];
		char* Query = argv[3];
		int numneighbors = atoi(argv[7]);
		string outfile = argv[4];
		int type = atoi(argv[5]);
		int wordlength = atoi(argv[6]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassTwoDistance(Alpha, Beta, Query, numneighbors, outfile, false, type, rand);
	}
	else if (argc == 10) {
		char* Tree = argv[1];
		char* BetaTree = argv[2];
		char* Query = argv[3];
		char type = *argv[4];
		double alphatau = atof(argv[5]);
		double betatau = atof(argv[6]);
		int numclasses = atoi(argv[7]);
		int wordlength = atoi(argv[8]);
		string outfile = argv[9];
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predict(Tree, BetaTree,  Query, type, alphatau, betatau , outfile, rand);
	}
	else if (argc == 11) {
		char* Tree = argv[1];
		char* BetaTree = argv[2];
		char* Query = argv[3];
		double alphatau = atof(argv[4]);
		double betatau = atof(argv[5]);
		int numclasses = atoi(argv[6]);
		int wordlength = atoi(argv[7]);
		string outfile = argv[8];
        int numneighbors = atoi(argv[9]);
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProtein(Tree, BetaTree,  Query, alphatau, betatau , outfile, numneighbors, rand);
	}
	else if (argc == 14) {
		char* Tree = argv[1];
		char* BetaTree = argv[2];
		char* Query = argv[3];
		double alphatau = atof(argv[4]);
		double betatau = atof(argv[5]);
		int numclasses = atoi(argv[6]);
		int wordlength = atoi(argv[7]);
		string outfile = argv[8];
        int numneighbors = atoi(argv[9]);
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProteinWithKNN(Tree, BetaTree,  Query, alphatau, betatau , outfile, numneighbors, rand);
	}
	/*else if (argc == 9){
		char* Tree = argv[1];
		char* Query = argv[2];
		int numneighbors = atoi(argv[3]);
		string outfile = argv[4];
		int numclasses = atoi(argv[5]);
		int wordlength = atoi(argv[6]);
		int firstIteration = atoi(argv[7]);
		bool isFirstIteration;
		if (firstIteration == 1)
			isFirstIteration = true;
		else
			isFirstIteration = false;
		int type = atoi(argv[8]);
		char* betatree = "../../code/filegen/beta.words";
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassDistance(Tree,betatree, Query, numneighbors, outfile, isFirstIteration, type);
	}
    */
	else if (argc == 5) {
		char* alpha = argv[1];
		char* beta = argv[2];
		char* alphaSaveTo = argv[3];
		char* betaSaveTo = argv[4];
		lpGenerator lpgen = lpGenerator(0,0,0,0,0,true);
		lpgen.buildSaveTrees(alpha, beta, alphaSaveTo, betaSaveTo);
	}
    else if (argc == 4){
        char* query = argv[1];
        char* neighbors = argv[2];
        string outfile = argv[3];
        lpGenerator lpgen = lpGenerator(0,0,0,0,0,true);
        lpgen.findDistances2(query,neighbors,outfile, rand);
    }
    else if (argc == 2){
        lpGenerator lpgen = lpGenerator(0,0,0,0,0,true);
        lpgen.printDFs();
    }

	else {
		ofstream output;
		string filename = "words.lp";
		output.open(filename.c_str());
		int numClasses = atoi( argv[1]);
		int k = atoi(argv[2]);
		int l = atoi(argv[3]);
		int TRSSize = atoi(argv[4]);
		int wordlength = atoi(argv[5]);
		char* alphats = argv[6];
		char* betats = argv[7];
		char* coilts = argv[8];

		lpGenerator lpgen = lpGenerator(numClasses, k, l, TRSSize, wordlength, true);

		lpgen.writeLP(alphats, betats, coilts, &output);
		cout << "I'm at the end of lpgen" << endl;
		output.close();
	}

}
