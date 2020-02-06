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
    //if (strcmp(argv[1],"FindNeighborsOneClass") == 0 || strcmp(argv[1],"CoverSet") == 0 || strcmp(argv[1],"CombinedNeighbors") == 0 || strcmp(argv[1],"OneClassTwoDistance") == 0 || strcmp(argv[1],"OneClassTwoDistancek") == 0 || strcmp(argv[1],"GreedyPredict") == 0 || strcmp(argv[1],"PredictProtein") == 0 || strcmp(argv[1],"PredictProteinKnn") == 0 || strcmp(argv[1],"FindDistances2") == 0){
        //ifstream src(argv[2],ios::binary);
        //ofstream dst("alphadis.h", ios::binary);
        //dst << src.rdbuf();
        //src.close();
        //dst.close();
        //ifstream src(argv[3],ios::binary);
        //ofstream dst("betadis.h", ios::binary);
        //dst << src.rdbuf();
        //src.close();
        //dst.close();
    //}

	if (strcmp(argv[1],"lpgen") == 0) {
		ofstream output;
		string filename = "words.lp";
		output.open(filename.c_str());
		int numClasses = atoi( argv[8]);
		int k = atoi(argv[9]);
		int l = atoi(argv[10]);
		int TRSSize = atoi(argv[11]);
		int wordlength = atoi(argv[12]);
		int firstIteration = atoi(argv[13]);
		bool isFirstIteration;
		if (firstIteration == 1)
			isFirstIteration = true;
		else
			isFirstIteration = false;
		lpGenerator lpgen = lpGenerator(numClasses, k, l, TRSSize, wordlength, false);
        rand = atoi(argv[2]);
		char* alpha = argv[3];
		char* beta = argv[4];
		char* alphats = argv[5];
		char* betats = argv[6];
		char* coilts = argv[7];
		cout << "starting to find neighbors" << endl;
		lpgen.findNeighbors(alpha, beta, alphats, betats, coilts, isFirstIteration, rand);
		cout << "found all the neighbors, now writing to file" << endl;
		lpgen.writeLP(alphats, betats, coilts, &output);
		cout << "I'm at the end of lpgen" << endl;
		output.close();
	}
    //From manyjobs.sh
    else if (strcmp(argv[1],"FindNeighborsOneClass") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* Query = argv[4];
		int numneighbors = atoi(argv[5]);
		string outfile = argv[6];
		int numclasses = atoi(argv[7]);
		int wordlength = atoi(argv[8]);
		int firstIteration = atoi(argv[9]);
		bool isFirstIteration;
        int type = atoi(argv[10]);
        float rho = atof(argv[11]);
		if (firstIteration == 1)
			isFirstIteration = true;
		else
			isFirstIteration = false;
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClass(Tree, Query, numneighbors, outfile, isFirstIteration,type, rand, rho, false, false);
	}
    else if (strcmp(argv[1],"FindNeighborsOneClassLoad") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* Query = argv[4];
		int numneighbors = atoi(argv[5]);
		string outfile = argv[6];
		int numclasses = atoi(argv[7]);
		int wordlength = atoi(argv[8]);
		int firstIteration = atoi(argv[9]);
		bool isFirstIteration;
        int type = atoi(argv[10]);
        float rho = atof(argv[11]);
		if (firstIteration == 1)
			isFirstIteration = true;
		else
			isFirstIteration = false;
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClass(Tree, Query, numneighbors, outfile, isFirstIteration,type, rand, rho, false, true);
	}
    else if (strcmp(argv[1],"FindNeighborsOneClassCounts") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* Query = argv[4];
		int numneighbors = atoi(argv[5]);
		string outfile = argv[6];
		int numclasses = atoi(argv[7]);
		int wordlength = atoi(argv[8]);
		int firstIteration = atoi(argv[9]);
		bool isFirstIteration;
        int type = atoi(argv[10]);
        float rho = atof(argv[11]);
        string Trainingfile = argv[12];
		if (firstIteration == 1)
			isFirstIteration = true;
		else
			isFirstIteration = false;
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
        lpgen.buildmap(Trainingfile);
		lpgen.findNeighborsOneClass(Tree, Query, numneighbors, outfile, isFirstIteration,type, rand, rho, true, false);
	}
    else if (strcmp(argv[1],"FindNeighborsOneClassCountsLoad") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* Query = argv[4];
		int numneighbors = atoi(argv[5]);
		string outfile = argv[6];
		int numclasses = atoi(argv[7]);
		int wordlength = atoi(argv[8]);
		int firstIteration = atoi(argv[9]);
		bool isFirstIteration;
        int type = atoi(argv[10]);
        float rho = atof(argv[11]);
        string Trainingfile = argv[12];
		if (firstIteration == 1)
			isFirstIteration = true;
		else
			isFirstIteration = false;
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
        lpgen.buildmap(Trainingfile);
		lpgen.findNeighborsOneClass(Tree, Query, numneighbors, outfile, isFirstIteration,type, rand, rho, true, true);
	}
    else if (strcmp(argv[1],"CombinedNeighbors") == 0) {
        cout << "I'm here in the find combined" << endl;
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* Query = argv[4];
		int numneighbors = atoi(argv[5]);
		string outfile = argv[6];
		int numclasses = atoi(argv[7]);
		int wordlength = atoi(argv[8]);
		int firstIteration = atoi(argv[9]);
		bool isFirstIteration;
        int type = atoi(argv[10]);
        float epsilon = atof(argv[11]);
        float aarho = atof(argv[12]);
        float ssrho = atof(argv[13]);
		if (firstIteration == 1)
			isFirstIteration = true;
		else
			isFirstIteration = false;
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassCombined(Tree, Query, numneighbors, outfile, isFirstIteration,type, rand, epsilon, aarho, ssrho);
	}
    else if (strcmp(argv[1],"CoverSet") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* Query = argv[4];
		int numneighbors = atoi(argv[5]);
		string outfile = argv[6];
		int numclasses = atoi(argv[7]);
		int wordlength = atoi(argv[8]);
		int firstIteration = atoi(argv[9]);
		bool isFirstIteration;
        int type = atoi(argv[10]);
        if (firstIteration == 1)
            isFirstIteration = true;
        else
            isFirstIteration = false;
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsCoverset(Tree, Query, numneighbors, outfile, isFirstIteration,type, rand);
	}
    //From twodist.sh
    else if (strcmp(argv[1],"OneClassTwoDistance") == 0) {
        cout << "argc was 7" << endl;
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
		char* Query = argv[5];
		int numneighbors = 2;
		string outfile = argv[6];
		int type = atoi(argv[7]);
		int wordlength = atoi(argv[8]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassTwoDistance(Alpha, Beta, Query, numneighbors, outfile, false, type, rand, false,false, false);
	}
    else if (strcmp(argv[1],"OneClassTwoDistanceAve") == 0) {
        cout << "argc was 7" << endl;
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
		char* Query = argv[5];
		int numneighbors = atoi(argv[9]) + 1;
		string outfile = argv[6];
		int type = atoi(argv[7]);
		int wordlength = atoi(argv[8]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassTwoDistanceAve(Alpha, Beta, Query, numneighbors, outfile, false, type, rand,false);
	}
    else if (strcmp(argv[1],"OneClassTwoDistanceAveLoad") == 0) {
        cout << "argc was 7" << endl;
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
		char* Query = argv[5];
		int numneighbors = atoi(argv[9]) + 1;
		string outfile = argv[6];
		int type = atoi(argv[7]);
		int wordlength = atoi(argv[8]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassTwoDistanceAveLoad(Alpha, Beta, Query, numneighbors, outfile, false, type, rand,false);
	}
    else if (strcmp(argv[1],"OneClassTwoDistanceAveLoadCounts") == 0) {
        cout << "argc was 7" << endl;
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
		char* Query = argv[5];
		int numneighbors = atoi(argv[9]);
		string outfile = argv[6];
		int type = atoi(argv[7]);
		int wordlength = atoi(argv[8]);
        string Trainingfile = argv[10];
		int firstIteration = atoi(argv[11]);
		bool isFirstIteration;
        if (firstIteration == 1)
            isFirstIteration = true;
        else
            isFirstIteration = false;

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
        lpgen.buildmap(Trainingfile);
		lpgen.findNeighborsOneClassTwoDistanceAveLoad(Alpha, Beta, Query, numneighbors, outfile, isFirstIteration, type, rand,true);
	}
    else if (strcmp(argv[1],"OneClassTwoDistanceAveLoadCountsCoil") == 0) {
        cout << "argc was 7" << endl;
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
        char* Coil = argv[5];
		char* Query = argv[6];
		int numneighbors = atoi(argv[10]);
		string outfile = argv[7];
		int type = atoi(argv[8]);
		int wordlength = atoi(argv[9]);
        string Trainingfile = argv[11];
		int firstIteration = atoi(argv[12]);
		bool isFirstIteration;
        if (firstIteration == 1)
            isFirstIteration = true;
        else
            isFirstIteration = false;

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
        lpgen.buildmap(Trainingfile);
		lpgen.findNeighborsOneClassTwoDistanceAveLoadCoil(Alpha, Beta,Coil, Query, numneighbors, outfile, isFirstIteration, type, rand,true);
	}
    else if (strcmp(argv[1],"OneClassTwoDistanceAveLoadCountsDefault") == 0) {
        cout << "argc was 7" << endl;
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
		char* Query = argv[5];
		int numneighbors = atoi(argv[9]);
		string outfile = argv[6];
		int type = atoi(argv[7]);
		int wordlength = atoi(argv[8]);
        string Trainingfile = argv[10];

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
        lpgen.buildmap(Trainingfile);
		lpgen.findNeighborsOneClassTwoDistanceAveLoad(Alpha, Beta, Query, numneighbors, outfile, true, type, rand,true);
	}
    else if (strcmp(argv[1],"OneClassTwoDistanceAveShift") == 0) {
        cout << "argc was 7" << endl;
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
		char* Query = argv[5];
		int numneighbors = atoi(argv[9]) + 1;
		string outfile = argv[6];
		int type = atoi(argv[7]);
		int wordlength = atoi(argv[8]);
        float shift = atof(argv[9]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
        lpgen.setshift(shift);
		lpgen.findNeighborsOneClassTwoDistanceAve(Alpha, Beta, Query, numneighbors, outfile, false, type, rand,false);
	}
    else if (strcmp(argv[1],"OneClassTwoDistanceWeighted") == 0) {
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
		char* Query = argv[5];
		int numneighbors = atoi(argv[9]);
		string outfile = argv[6];
		int type = atoi(argv[7]);
		int wordlength = atoi(argv[8]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassTwoDistance(Alpha, Beta, Query, numneighbors, outfile, false, type, rand, false, true,false);
	}
    else if (strcmp(argv[1],"OneClassTwoDistanceWeightedLoad") == 0) {
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
		char* Query = argv[5];
		int numneighbors = atoi(argv[9]);
		string outfile = argv[6];
		int type = atoi(argv[7]);
		int wordlength = atoi(argv[8]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassTwoDistance(Alpha, Beta, Query, numneighbors, outfile, false, type, rand, true, true, false);
	}
    else if (strcmp(argv[1],"OneClassTwoDistanceWeightedDefaultLoad") == 0) {
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
		char* Query = argv[5];
		int numneighbors = atoi(argv[9]);
		string outfile = argv[6];
		int type = atoi(argv[7]);
		int wordlength = atoi(argv[8]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassTwoDistance(Alpha, Beta, Query, numneighbors, outfile, true, type, rand, true, true, false);
	}
    else if (strcmp(argv[1],"OneClassTwoDistanceWeightedDefaultLoadCoil") == 0) {
        rand = atoi(argv[2]);
		char* Coil = argv[3];
		char* Query = argv[4];
		int numneighbors = atoi(argv[8]);
		string outfile = argv[5];
		int type = atoi(argv[6]);
		int wordlength = atoi(argv[7]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsCoilClassTwoDistance(Coil, Query, numneighbors, outfile, true, type, rand, true, true, false);
	}
    else if (strcmp(argv[1],"OneClassTwoDistancek") == 0) {
        cout << "argc was 8" << endl;
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
		char* Query = argv[5];
		int numneighbors = atoi(argv[9]);
		string outfile = argv[6];
		int type = atoi(argv[7]);
		int wordlength = atoi(argv[8]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassTwoDistance(Alpha, Beta, Query, numneighbors, outfile, false, type, rand, false,false,false);
	}
    else if (strcmp(argv[1],"OneClassTwoDistancekDefault") == 0) {
        cout << "argc was 8" << endl;
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
		char* Query = argv[5];
		int numneighbors = atoi(argv[9]);
		string outfile = argv[6];
		int type = atoi(argv[7]);
		int wordlength = atoi(argv[8]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassTwoDistance(Alpha, Beta, Query, numneighbors, outfile, true, type, rand, false,false, false);
	}
    else if (strcmp(argv[1],"OneClassTwoDistancekDefaultLoad") == 0) {
        cout << "argc was 8" << endl;
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
		char* Query = argv[5];
		int numneighbors = atoi(argv[9]);
		string outfile = argv[6];
		int type = atoi(argv[7]);
		int wordlength = atoi(argv[8]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassTwoDistance(Alpha, Beta, Query, numneighbors, outfile, true, type, rand, true,false,false);
	}
    else if (strcmp(argv[1],"OneClassTwoDistanceAveDefault") == 0) {
        cout << "argc was 8" << endl;
        rand = atoi(argv[2]);
		char* Alpha = argv[3];
		char* Beta = argv[4];
		char* Query = argv[5];
		int numneighbors = atoi(argv[9]);
		string outfile = argv[6];
		int type = atoi(argv[7]);
		int wordlength = atoi(argv[8]);

		lpGenerator lpgen = lpGenerator(3, 0, 0, 0, wordlength, true);
		lpgen.findNeighborsOneClassTwoDistanceAve(Alpha, Beta, Query, numneighbors, outfile, true, type, rand,false);
	}
    else if (strcmp(argv[1],"GreedyPredict") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* BetaTree = argv[4];
		char* Query = argv[5];
		char type = *argv[6];
		double alphatau = atof(argv[7]);
		double betatau = atof(argv[8]);
		int numclasses = atoi(argv[9]);
		int wordlength = atoi(argv[10]);
		string outfile = argv[11];
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predict(Tree, BetaTree,  Query, type, alphatau, betatau , outfile, rand);
	}
    else if (strcmp(argv[1],"PredictProtein") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* BetaTree = argv[4];
		char* Query = argv[5];
		double alphatau = atof(argv[6]);
		double betatau = atof(argv[7]);
		int numclasses = atoi(argv[8]);
		int wordlength = atoi(argv[9]);
		string outfile = argv[10];
        int numneighbors = atoi(argv[11]);
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProtein(Tree, BetaTree,  Query, alphatau, betatau , outfile, numneighbors, rand,false);
	}
    else if (strcmp(argv[1],"PredictProteinDefault") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* BetaTree = argv[4];
		char* Query = argv[5];
		double alphatau = atof(argv[6]);
		double betatau = atof(argv[7]);
		int numclasses = atoi(argv[8]);
		int wordlength = atoi(argv[9]);
		string outfile = argv[10];
        int numneighbors = atoi(argv[11]);
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProtein(Tree, BetaTree,  Query, alphatau, betatau , outfile, numneighbors, rand,true);
	}
    else if (strcmp(argv[1],"PredictProteinLoadTree") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* BetaTree = argv[4];
		char* Query = argv[5];
		double alphatau = atof(argv[6]);
		double betatau = atof(argv[7]);
		int numclasses = atoi(argv[8]);
		int wordlength = atoi(argv[9]);
		string outfile = argv[10];
        int numneighbors = atoi(argv[11]);
        string treefile = "trees/alpha" + to_string(rand) + ".tree";
        string treefile2 = "trees/beta" + to_string(rand) + ".tree";
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProteinLoadTree(Tree, BetaTree,  Query, alphatau, betatau , outfile, numneighbors, rand, treefile, treefile2, false, false);
	}
    else if (strcmp(argv[1],"PredictProteinLoadTreeWeights") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* BetaTree = argv[4];
		char* Query = argv[5];
		double alphatau = atof(argv[6]);
		double betatau = atof(argv[7]);
		int numclasses = atoi(argv[8]);
		int wordlength = atoi(argv[9]);
		string outfile = argv[10];
        int numneighbors = atoi(argv[11]);
        string treefile = "trees/alpha" + to_string(rand) + ".tree";
        string treefile2 = "trees/beta" + to_string(rand) + ".tree";
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProteinLoadTree(Tree, BetaTree,  Query, alphatau, betatau , outfile, numneighbors, rand, treefile, treefile2, false, true);
	}
    else if (strcmp(argv[1],"PredictProteinLoadTreeDefault") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* BetaTree = argv[4];
		char* Query = argv[5];
		double alphatau = atof(argv[6]);
		double betatau = atof(argv[7]);
		int numclasses = atoi(argv[8]);
		int wordlength = atoi(argv[9]);
		string outfile = argv[10];
        int numneighbors = atoi(argv[11]);
        string treefile = "trees/alpha" + to_string(rand) + argv[12] + ".tree";
        string treefile2 = "trees/beta" + to_string(rand) + argv[12] + ".tree";
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProteinLoadTreeDefault(Tree, BetaTree,  Query, alphatau, betatau , outfile, numneighbors, rand, treefile,treefile2, atoi(argv[13]), false, false);
	}
    else if (strcmp(argv[1],"PredictProteinLoadTreeDefaultWeights") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* BetaTree = argv[4];
		char* Query = argv[5];
		double alphatau = atof(argv[6]);
		double betatau = atof(argv[7]);
		int numclasses = atoi(argv[8]);
		int wordlength = atoi(argv[9]);
		string outfile = argv[10];
        int numneighbors = atoi(argv[11]);
        string treefile = "trees/alpha" + to_string(rand) + argv[12] + ".tree";
        string treefile2 = "trees/beta" + to_string(rand) + argv[12] + ".tree";
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProteinLoadTreeDefault(Tree, BetaTree,  Query, alphatau, betatau , outfile, numneighbors, rand, treefile,treefile2, atoi(argv[13]), true, false);
	}
    else if (strcmp(argv[1],"PredictProteinLoadTreeDefaultWeightsCoil") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* BetaTree = argv[4];
		char* CoilTree = argv[5];
		char* Query = argv[6];
		double alphatau = atof(argv[7]);
		double betatau = atof(argv[8]);
		int numclasses = atoi(argv[9]);
		int wordlength = atoi(argv[10]);
		string outfile = argv[11];
        int numneighbors = atoi(argv[12]);
        string treefile = "trees/alpha" + to_string(rand) + argv[13] + ".tree";
        string treefile2 = "trees/beta" + to_string(rand) + argv[13] + ".tree";
        string treefile3 = "trees/coil" + to_string(rand) + argv[12] + ".tree";
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProteinLoadTreeDefaultCoil(Tree, BetaTree,  CoilTree,Query, alphatau, betatau , outfile, numneighbors, rand, treefile,treefile2,treefile3, atoi(argv[13]), true, false);
    }
    else if (strcmp(argv[1],"PredictProteinLoadTreeCoil") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* BetaTree = argv[4];
		char* CoilTree = argv[5];
		char* Query = argv[6];
		double alphatau = atof(argv[7]);
		double betatau = atof(argv[8]);
		int numclasses = atoi(argv[9]);
		int wordlength = atoi(argv[10]);
		string outfile = argv[11];
        int numneighbors = atoi(argv[12]);
        string treefile = "trees/alpha" + to_string(rand) + ".tree";
        string treefile2 = "trees/beta" + to_string(rand) + ".tree";
        string treefile3 = "trees/coil" + to_string(rand) + ".tree";
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProteinLoadTreeCoil(Tree, BetaTree,  CoilTree,Query, alphatau, betatau , outfile, numneighbors, rand, treefile,treefile2,treefile3, atoi(argv[13]), false, false);
	}
    else if (strcmp(argv[1],"TrainingLoadTree") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* BetaTree = argv[4];
		char* Query = argv[5];
		double alphatau = atof(argv[6]);
		double betatau = atof(argv[7]);
		int numclasses = atoi(argv[8]);
		int wordlength = atoi(argv[9]);
		string outfile = argv[10];
        int numneighbors = atoi(argv[11]);
        string treefile = "trees/alpha" + to_string(rand) + argv[12] + ".tree";
        string treefile2 = "trees/beta" + to_string(rand) + argv[12] + ".tree";
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProteinLoadTreeDefault(Tree, BetaTree,  Query, alphatau, betatau , outfile, numneighbors, rand, treefile,treefile2,atoi(argv[13]), true, true);
	}
    else if (strcmp(argv[1],"PredictProteinKnn") == 0) {
        rand = atoi(argv[2]);
		char* Tree = argv[3];
		char* BetaTree = argv[4];
		char* Query = argv[5];
		double alphatau = atof(argv[6]);
		double betatau = atof(argv[7]);
		int numclasses = atoi(argv[8]);
		int wordlength = atoi(argv[9]);
		string outfile = argv[10];
        int numneighbors = atoi(argv[11]);
		lpGenerator lpgen = lpGenerator(numclasses, 0, 0, 0, wordlength, true);
		lpgen.predictProteinWithKNN(Tree, BetaTree,  Query, alphatau, betatau , outfile, numneighbors, rand);
	}
    else if (strcmp(argv[1],"PredictProteinKnnCoil") == 0) {
        rand = atoi(argv[2]);
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
		lpgen.predictProteinWithKNNCoil(Tree, BetaTree,CoilTree,  Query , outfile, numneighbors, rand,isFirstIteration, dirnum,false,false,false);
    }
    else if (strcmp(argv[1],"PredictProteinKnnCoilDirr") == 0) {
        printf("predictproteinknncoildirr");
        rand = atoi(argv[2]);
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
		lpgen.predictProteinWithKNNCoil(Tree, BetaTree,CoilTree,  Query , outfile, numneighbors, rand,isFirstIteration, dirnum,false,false,true);
    }
    else if (strcmp(argv[1],"PredictProteinKnnCoilWeights") == 0) {
        rand = atoi(argv[2]);
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
		lpgen.predictProteinWithKNNCoil(Tree, BetaTree,CoilTree,  Query , outfile, numneighbors, rand,isFirstIteration, dirnum,false,true,false);
    }
    else if (strcmp(argv[1],"SSPredictProteinKnnCoil") == 0) {
        rand = atoi(argv[2]);
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
		lpgen.predictProteinWithKNNCoil(Tree, BetaTree,CoilTree,  Query , outfile, numneighbors, rand,isFirstIteration, dirnum,true,false,false);
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
    else if (strcmp(argv[1],"SaveTrees") == 0) {
        int rand = atoi(argv[2]);
		char* alpha = argv[3];
		char* beta = argv[4];
		char* alphaSaveTo = argv[5];
		char* betaSaveTo = argv[6];
		lpGenerator lpgen = lpGenerator(0,0,0,0,23,true);
		lpgen.buildSaveTrees(alpha, beta, alphaSaveTo, betaSaveTo, rand,false);
	}
    else if (strcmp(argv[1],"SaveTreesCoil") == 0) {
        int rand = atoi(argv[2]);
		char* alpha = argv[3];
		char* beta = argv[4];
        char* coil = argv[5];
		char* alphaSaveTo = argv[6];
		char* betaSaveTo = argv[7];
		char* coilSaveTo = argv[8];
		lpGenerator lpgen = lpGenerator(0,0,0,0,23,true);
		lpgen.buildSaveTreesCoil(alpha, beta,coil, alphaSaveTo, betaSaveTo,coilSaveTo, rand,false);
	}
    else if (strcmp(argv[1],"SaveTreesCoilWeights") == 0) {
        int rand = atoi(argv[2]);
		char* alpha = argv[3];
		char* beta = argv[4];
        char* coil = argv[5];
		char* alphaSaveTo = argv[6];
		char* betaSaveTo = argv[7];
		char* coilSaveTo = argv[8];
		lpGenerator lpgen = lpGenerator(0,0,0,0,23,true);
		lpgen.buildSaveTreesCoil(alpha, beta,coil, alphaSaveTo, betaSaveTo,coilSaveTo, rand,true);
	}
    else if (strcmp(argv[1],"SaveSSTrees") == 0) {
        int rand = atoi(argv[2]);
		char* alpha = argv[3];
		char* beta = argv[4];
        char* coil = argv[5];
		char* alphaSaveTo = argv[6];
		char* betaSaveTo = argv[7];
		char* coilSaveTo = argv[8];
		lpGenerator lpgen = lpGenerator(0,0,0,0,23,true);
		lpgen.buildSStrees(alpha, beta,coil, alphaSaveTo, betaSaveTo,coilSaveTo, rand,true);
        //TODO: fix so that if there aren't any weights we don't make a weighted tree
	}
    else if (strcmp(argv[1],"SaveTreesWeights") == 0) {
        int rand = atoi(argv[2]);
		char* alpha = argv[3];
		char* beta = argv[4];
		char* alphaSaveTo = argv[5];
		char* betaSaveTo = argv[6];
		lpGenerator lpgen = lpGenerator(0,0,0,0,23,true);
		lpgen.buildSaveTrees(alpha, beta, alphaSaveTo, betaSaveTo, rand,true);
	}
    else if (strcmp(argv[1],"SaveTreesDefault") == 0) {
        int rand = atoi(argv[2]);
		char* alpha = argv[3];
		char* beta = argv[4];
		char* alphaSaveTo = argv[5];
		char* betaSaveTo = argv[6];
		lpGenerator lpgen = lpGenerator(0,0,0,0,23,true);
		lpgen.buildSaveTreesDefault(alpha, beta, alphaSaveTo, betaSaveTo, rand,false);
	}
    else if (strcmp(argv[1],"SaveTreesDefaultWeighted") == 0) {
        int rand = atoi(argv[2]);
		char* alpha = argv[3];
		char* beta = argv[4];
		char* alphaSaveTo = argv[5];
		char* betaSaveTo = argv[6];
		lpGenerator lpgen = lpGenerator(0,0,0,0,23,true);
		lpgen.buildSaveTreesDefault(alpha, beta, alphaSaveTo, betaSaveTo, rand,true);
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
    else if (strcmp(argv[1],"FindDistances2") == 0) {
        rand = atoi(argv[2]);
        char* query = argv[3];
        char* neighbors = argv[4];
        string outfile = argv[5];
        lpGenerator lpgen = lpGenerator(0,0,0,0,23,true);
        lpgen.findDistances2(query,neighbors,outfile, rand);
    }
    else if (strcmp(argv[1],"FindCombinedDistances2") == 0) {
        rand = atoi(argv[2]);
        char* query = argv[3];
        char* neighbors = argv[4];
        string outfile = argv[5];
        float epsilon = atof(argv[6]);
        float aarho = atof(argv[7]);
        float ssrho = atof(argv[8]);
        lpGenerator lpgen = lpGenerator(0,0,0,0,23,true);
        lpgen.findCombinedDistances2(query,neighbors,outfile, rand, epsilon,aarho, ssrho);
    }
    else if (strcmp(argv[1],"PrintDFs") == 0) {
        lpGenerator lpgen = lpGenerator(0,0,0,0,0,true);
        lpgen.printDFs();
    }

    else if (strcmp(argv[1],"WriteLP") == 0) {
		ofstream output;
		string filename = "words.lp";
		output.open(filename.c_str());
		int numClasses = atoi( argv[2]);
		int k = atoi(argv[3]);
		int l = atoi(argv[4]);
		int TRSSize = atoi(argv[5]);
		int wordlength = atoi(argv[6]);
		char* alphats = argv[7];
		char* betats = argv[8];
		char* coilts = argv[9];

		lpGenerator lpgen = lpGenerator(numClasses, k, l, TRSSize, wordlength, true);

		lpgen.writeLP(alphats, betats, coilts, &output);
		cout << "I'm at the end of lpgen" << endl;
		output.close();
	}
    else if (strcmp(argv[1],"ct") == 0) {
        lpGenerator lpgen = lpGenerator(0,0,0,0,23,true);
        lpgen.gettargetsdist(5);
    }

    else{
        lpGenerator lpgen = lpGenerator(3,3,3,3,23,true);
        lpgen.getnewtargets(atoi(argv[1]));
        //lpgen.testDF(atoi(argv[1]),argv[2],argv[3]);
        cout << "Error in function call" << endl;
    }
    return 0;

}
