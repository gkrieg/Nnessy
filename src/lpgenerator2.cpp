/**@file	/gsfs1/xdisk/skrieger/code/mgenerator/lpgenerator.cpp
 * @author	skrieger
 * @version	702
 * @date
 * 	Created:	Thu 23 Jun 2016 02:10:10 PM MST \n
 * 	Last Update:	Tue Jan 29 2019 02:10:10 PM MST
 */

/*===========================================================================*/
/*===============================[ ]===============================*/
/*===========================================================================*/
#include <string>
#include "lpgenerator.h"
#include <fstream>
#include <iostream>
#include <map>
#include <unistd.h>
#include <thread>
#include <stdlib.h>
#include <string.h>
using namespace std;
void lpGenerator::buildSaveTreesCoilDefault(char* alpha, char* beta,char* coil, char* alphaSave, char* betaSave,char* coilSave,int rand,bool weights) {
  const char* AlphaInputFile = (const char* ) alpha;
  const char* BetaInputFile = (const char* ) beta;
  const char* CoilInputFile = (const char* ) coil;
  string r = "r";
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  FILE* coilFile = fopen(CoilInputFile, r.c_str());
  setEpsilonRho(1,0.9,0.3);
  if (weights){
      string weightsfilename = "weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
      thread alphatree(&lpGenerator::buildalphatreeweighted, this, alphaFile, 1, true, rand);
      thread betatree(&lpGenerator::buildbetatreeweighted, this, betaFile, 1, true, rand);
      thread coiltree(&lpGenerator::buildcoiltreeweighted, this, coilFile, 1, true, rand);
      alphatree.join();
      betatree.join();
      coiltree.join();
  }else{
      thread alphatree(&lpGenerator::buildalphatree, this, alphaFile, 1, true, rand);
      thread betatree(&lpGenerator::buildbetatree, this, betaFile, 1, true, rand);
      thread coiltree(&lpGenerator::buildcoiltree, this, coilFile, 1, true, rand);
      alphatree.join();
      betatree.join();
      coiltree.join();
  }
  nnfinder1.saveDispersionTree(alphaSave, "alpha");
  nnfinder4.saveDispersionTree(betaSave, "beta");
  nnfinder7.saveDispersionTree(coilSave, "coil");
}
void lpGenerator::predictProteinWithKNNCoil(char* AlphaIF, char* BetaIF,char* CoilIF, char* ProteinIF, string filename, int numneighbors, int rand, bool firstIteration, int dirnum) {

  const char* AlphaInputFile = (const char* ) AlphaIF;
  const char* BetaInputFile = (const char* ) BetaIF;
  const char* CoilInputFile = (const char* ) CoilIF;
  string r = "r";
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  FILE* coilFile = fopen(CoilInputFile, r.c_str());
  cout << "###################################" << BetaIF << endl;
  Distance alphadistfunc;
  Distance betadistfunc;
  Distance coildistfunc;
  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  FILE* ProteinTSFile = fopen(ProteinIF, r.c_str());
  //at this point is when we need to split up the file and make each word into its own thing.  This is going to be the hardest part.
  MetricSpace* MSAlpha = CreateMetricSpace(inputchar, ProteinTSFile);
  lpGenerator::TRSSize = MSAlpha->numpoints;

  if (firstIteration){
      string weightsfilename = "weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
      alphadistfunc = WeightedWordDist;
      betadistfunc = WeightedWordDist;
      coildistfunc = WeightedWordDist;
      nnfinder1.loadDispersionTree(alphaFile, "trees/alpha" + to_string(rand) + to_string(dirnum) + ".tree", "defaultaw",rand);
      nnfinder1.queryAlpha(numneighbors, MSAlpha, &this->AANeighbors, 1);
      nnfinder1.destroyDTree('A');
      nnfinder4.loadDispersionTree(betaFile, "trees/beta" + to_string(rand) + to_string(dirnum) + ".tree","defaultbw",rand);
      nnfinder4.queryBeta(numneighbors, MSAlpha, &this->BBNeighbors, 1);
      nnfinder4.destroyDTree('B');
      printf("Destroyed the beta tree\n");
      nnfinder7.loadDispersionTree(coilFile, "trees/coil" + to_string(rand) + to_string(dirnum) + ".tree","defaultcw",rand);
      nnfinder7.queryCoil(numneighbors, MSAlpha, &this->CCNeighbors, 1);
      printf("before Destroyed the coil tree\n");
      nnfinder7.destroyDTree('C');
      printf("Destroyed the coil tree\n");
  }
  else{
      nnfinder1.loadDispersionTree(alphaFile, "trees/alpha" + to_string(rand) + ".tree", "alpha",rand);
      nnfinder4.loadDispersionTree(betaFile, "trees/beta" + to_string(rand) + ".tree","beta",rand);
      nnfinder7.loadDispersionTree(coilFile, "trees/coil" + to_string(rand) + ".tree","coil",rand);
      alphadistfunc = wordDist2Alpha;
      betadistfunc = wordDist2Beta;
      coildistfunc = wordDist2Coil;
      nnfinder1.queryAlpha(numneighbors, MSAlpha, &this->AANeighbors, 1);
      nnfinder4.queryBeta(numneighbors, MSAlpha, &this->BBNeighbors, 1);
      nnfinder7.queryCoil(numneighbors, MSAlpha, &this->CCNeighbors, 1);
  }

  //int numneighbors = 2;
  ofstream outfile(filename.c_str());
  float alphadist;
  float betadist;
  float coildist;
  for (int i = 0; i < TRSSize; i++) {
      for (int j = 0;j < numneighbors; j++){
		   alphadist = alphadistfunc(MSAlpha->points[i], AANeighbors[i][j]);
		   betadist = betadistfunc(MSAlpha->points[i], BBNeighbors[i][j]);
		   coildist = coildistfunc(MSAlpha->points[i], CCNeighbors[i][j]);
           outfile << (char*)MSAlpha->points[i] << " " << (char*)AANeighbors[i][j] << " " << alphadist << " " << (char*)BBNeighbors[i][j] << " " << betadist << " " << (char*)CCNeighbors[i][j] << " " << coildist << endl;
      }

	}

	outfile.close();
}



