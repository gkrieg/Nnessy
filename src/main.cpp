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
    if (strcmp(argv[1],"PredictProteinKnnCoil") == 0) {
        int rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* BetaTree = argv[4];
		char* CoilTree = argv[5];
		char* Query = argv[6];
        int numclasses = 3;
		int wordlength = atoi(argv[7]);
		string outfile = argv[8];
        int numneighbors = atoi(argv[9]);
		int firstIteration = atoi(argv[10]);
        int dirnum = atoi(argv[11]);
		bool isFirstIteration;
        if (firstIteration == 1)
            isFirstIteration = true;
        else
            isFirstIteration = false;
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProteinWithKNNCoil(Tree, BetaTree,CoilTree,  Query , outfile, numneighbors, rand,isFirstIteration, dirnum);
    }
    else if (strcmp(argv[1],"SaveTreesDefaultCoilWeighted") == 0) {
        int rand = atoi(argv[2]);
		char* alpha = argv[3];
		char* beta = argv[4];
		char* coil = argv[5];
		char* alphaSaveTo = argv[6];
		char* betaSaveTo = argv[7];
		char* coilSaveTo = argv[8];
        int wordlen = atoi(argv[9]);
		lpGenerator lpgen = lpGenerator(0,0,0,0,wordlen,true);
		lpgen.buildSaveTreesCoilDefault(alpha, beta,coil, alphaSaveTo, betaSaveTo,coilSaveTo, rand,true);
	}
    else{
        cout << "Error in function call" << endl;
    }
    return 0;

}
