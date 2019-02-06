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
		char* Tree = argv[2];
		char* BetaTree = argv[3];
		char* CoilTree = argv[4];
		char* Query = argv[5];
        int numclasses = 3;
		int wordlength = atoi(argv[6]);
		string outfile = argv[7];
        int numneighbors = atoi(argv[8]);
		int firstIteration = atoi(argv[9]);
		bool isFirstIteration;
        if (firstIteration == 1)
            isFirstIteration = true;
        else
            isFirstIteration = false;
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProteinWithKNNCoil(Tree, BetaTree,CoilTree,  Query , outfile, numneighbors,isFirstIteration);
    }
    else{
        cout << "Error in function call" << endl;
    }
    return 0;

}
