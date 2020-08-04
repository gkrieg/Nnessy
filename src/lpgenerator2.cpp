/**@file	/gsfs1/xdisk/skrieger/code/mgenerator/lpgenerator.cpp
 * @author	skrieger
 * @version	702
 * @date
 * 	Created:	Thu 23 Jun 2016 02:10:10 PM MST \n
 * 	Last Update:	Thu 23 Jun 2016 02:10:10 PM MST
 */

/*===========================================================================*/
/*===============================[ ]===============================*/
/*===========================================================================*/
#include <string>
#include "lpgenerator.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <thread>
#include <stdlib.h>
#include <string.h>
using namespace std;
const char *residues2 = "ARNDCQEGHILKMFPSTWYVBZX-";


void lpGenerator::predictProteinWithKNNCoil(char* AlphaIF, char* BetaIF,char* CoilIF, char* ProteinIF, string filename, int numneighbors, int rand, bool firstIteration, int dirnum,bool ss,bool weights,bool dirr) {

  const char* AlphaInputFile = (const char* ) AlphaIF;
  const char* BetaInputFile = (const char* ) BetaIF;
  const char* CoilInputFile = (const char* ) CoilIF;
  string r = "r";
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  FILE* coilFile = fopen(CoilInputFile, r.c_str());
  Distance alphadistfunc;
  Distance betadistfunc;
  Distance coildistfunc;
  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';
  int alphann = 0;
  int betann = 0;
  int coilnn = 0;

  FILE* ProteinTSFile = fopen(ProteinIF, r.c_str());
  //at this point is when we need to split up the file and make each word into its own thing.  This is going to be the hardest part.
  MetricSpace* MSAlpha = CreateMetricSpace(inputchar, ProteinTSFile);
  lpGenerator::TRSSize = MSAlpha->numpoints;

  if (firstIteration){
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
      if (ss){
          alphadistfunc = SecondaryDist2;
          betadistfunc = SecondaryDist2;
          coildistfunc = SecondaryDist2;
          nnfinder1.loadDispersionTree(alphaFile, "trees/ssalpha" + to_string(rand)  + ".tree", "ass",rand);
          nnfinder4.loadDispersionTree(betaFile, "trees/ssbeta" + to_string(rand) + ".tree","bss",rand);
          nnfinder7.loadDispersionTree(coilFile, "trees/sscoil" + to_string(rand) + ".tree","css",rand);
          int treesize = 12000000;
          alphann = nnfinder1.queryAlphaRange(treesize, MSAlpha, &this->AANeighbors, 1);
          nnfinder1.destroyDTree('A');
          betann = nnfinder4.queryBetaRange(treesize, MSAlpha, &this->BBNeighbors, 1);
          nnfinder4.destroyDTree('B');
          coilnn = nnfinder7.queryCoilRange(treesize, MSAlpha, &this->CCNeighbors, 1);
          nnfinder7.destroyDTree('C');
         /* 
          nnfinder1.loadDispersionTree(alphaFile, "a.tree", "ass",rand);
          nnfinder4.loadDispersionTree(betaFile, "b.tree","bss",rand);
          nnfinder7.loadDispersionTree(coilFile, "c.tree","css",rand);
          */
      } else {
          alphadistfunc = WeightedWordDist2;
          betadistfunc = WeightedWordDist2;
          coildistfunc = WeightedWordDist2;
          nnfinder1.loadDispersionTree(alphaFile, "trees/alpha.tree", "defaultaw",rand);
          nnfinder4.loadDispersionTree(betaFile, "trees/beta.tree","defaultbw",rand);
          nnfinder7.loadDispersionTree(coilFile, "trees/coil.tree","defaultcw",rand);
          nnfinder1.queryAlpha(numneighbors, MSAlpha, &this->AANeighbors, 1);
          nnfinder1.destroyDTree('A');
          nnfinder4.queryBeta(numneighbors, MSAlpha, &this->BBNeighbors, 1);
          nnfinder4.destroyDTree('B');
          nnfinder7.queryCoil(numneighbors, MSAlpha, &this->CCNeighbors, 1);
          nnfinder7.destroyDTree('C');
      }
  }
  else{
      if(weights){
          string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
          readinweights(this->lengthOfWord,weightsfilename.c_str());
          nnfinder1.loadDispersionTree(alphaFile, "trees/walpha" + to_string(rand) + ".tree", "alphaweighted",rand);
          nnfinder4.loadDispersionTree(betaFile, "trees/wbeta" + to_string(rand) + ".tree","betaweighted",rand);
          nnfinder7.loadDispersionTree(coilFile, "trees/wcoil" + to_string(rand) + ".tree","coilweighted",rand);
          alphadistfunc = WeightedAlphaDistance;
          betadistfunc = WeightedBetaDistance;
          coildistfunc = WeightedCoilDistance;
      }
      else{
          if(dirr){
          nnfinder1.loadDispersionTree(alphaFile, "trees/alpha" + to_string(dirnum) + to_string(rand) + ".tree", "alpha",rand);
          nnfinder4.loadDispersionTree(betaFile, "trees/beta" + to_string(dirnum) + to_string(rand) + ".tree","beta",rand);
          nnfinder7.loadDispersionTree(coilFile, "trees/coil" + to_string(dirnum) + to_string(rand) + ".tree","coil",rand);
          alphadistfunc = wordDist2Alpha;
          betadistfunc = wordDist2Beta;
          coildistfunc = wordDist2Coil;
          }
          else{
          nnfinder1.loadDispersionTree(alphaFile, "trees/alpha.tree", "alpha",rand);
          nnfinder4.loadDispersionTree(betaFile, "trees/beta.tree","beta",rand);
          nnfinder7.loadDispersionTree(coilFile, "trees/coil.tree","coil",rand);
          alphadistfunc = wordDist2Alpha;
          betadistfunc = wordDist2Beta;
          coildistfunc = wordDist2Coil;
          }
      }
      nnfinder1.queryAlpha(numneighbors, MSAlpha, &this->AANeighbors, 1);
      nnfinder4.queryBeta(numneighbors, MSAlpha, &this->BBNeighbors, 1);
      nnfinder7.queryCoil(numneighbors, MSAlpha, &this->CCNeighbors, 1);
      nnfinder1.destroyDTree('A');
      nnfinder4.destroyDTree('B');
      nnfinder7.destroyDTree('C');
  }

  //int numneighbors = 2;
  ofstream outfile(filename.c_str());
  float alphadist;
  float betadist;
  float coildist;
  float salphadist;
  float sbetadist;
  float scoildist;
  Distance adistfunc;
  Distance bdistfunc;
  Distance cdistfunc;
  if (ss){
      adistfunc = SecondaryDist2;
      bdistfunc = SecondaryDist2;
      cdistfunc = SecondaryDist2;
      alphadistfunc = WeightedWordDist2;
      betadistfunc = WeightedWordDist2;
      coildistfunc = WeightedWordDist2;
      for (int i = 0; i < TRSSize; i++) {
          vector<tuple<float,Pointer,float>> v;
          int j = 0;
          while(strcmp((char*)AANeighbors[i][j],"null") != 0){

               alphadist = alphadistfunc(MSAlpha->points[i], AANeighbors[i][j]);
               salphadist = adistfunc(MSAlpha->points[i], AANeighbors[i][j]);
               v.push_back(make_tuple(alphadist,AANeighbors[i][j],salphadist));

               //outfile << (char*)MSAlpha->points[i] << " " << (char*)AANeighbors[i][j] << " " << alphadist << " " <<salphadist <<endl;
              j++;
          }
          j = 0;
          while(strcmp((char*)BBNeighbors[i][j],"null") != 0){
               betadist = betadistfunc(MSAlpha->points[i], BBNeighbors[i][j]);
               sbetadist = bdistfunc(MSAlpha->points[i], BBNeighbors[i][j]);
               v.push_back(make_tuple(betadist,BBNeighbors[i][j],sbetadist));
               //outfile << (char*)MSAlpha->points[i] << " " << (char*)BBNeighbors[i][j] << " " << betadist << " " << sbetadist << endl;
              j++;
          }
          j = 0;
          while(strcmp((char*)CCNeighbors[i][j],"null") != 0){
               coildist = coildistfunc(MSAlpha->points[i], CCNeighbors[i][j]);
               scoildist = cdistfunc(MSAlpha->points[i], CCNeighbors[i][j]);
               v.push_back(make_tuple(coildist,CCNeighbors[i][j],scoildist));
               //outfile << (char*)MSAlpha->points[i] << " "<< (char*)CCNeighbors[i][j] << " " << coildist << " " << scoildist << endl;
              j++;
          }
          sort(v.begin(),v.end());
          int vectorend = v.size() > numneighbors ? numneighbors : v.size();
          if (vectorend == 0)
              outfile << (char*)MSAlpha->points[i] << endl;
          for (int r = 0;r < vectorend;r++){
              alphadist = get<0>(v[r]);
              salphadist = get<2>(v[r]);
               outfile << (char*)MSAlpha->points[i] << " " << (char*)get<1>(v[r]) << " " << alphadist << " " <<salphadist <<endl;
      }
      }


    }
  else{
  for (int i = 0; i < TRSSize; i++) {
      for (int j = 0;j < numneighbors; j++){
		   alphadist = alphadistfunc(MSAlpha->points[i], AANeighbors[i][j]);
		   betadist = betadistfunc(MSAlpha->points[i], BBNeighbors[i][j]);
		   coildist = coildistfunc(MSAlpha->points[i], CCNeighbors[i][j]);
           if (ss){
               salphadist = adistfunc(MSAlpha->points[i], AANeighbors[i][j]);
               sbetadist = bdistfunc(MSAlpha->points[i], BBNeighbors[i][j]);
               scoildist = cdistfunc(MSAlpha->points[i], CCNeighbors[i][j]);

               outfile << (char*)MSAlpha->points[i] << " " << (char*)AANeighbors[i][j] << " " << alphadist << " " <<salphadist << " " << (char*)BBNeighbors[i][j] << " " << betadist << " " << sbetadist << " "<< (char*)CCNeighbors[i][j] << " " << coildist << " " << scoildist << endl;
           }else{
               outfile << (char*)MSAlpha->points[i] << " " << (char*)AANeighbors[i][j] << " " << alphadist << " " << (char*)BBNeighbors[i][j] << " " << betadist << " "<< (char*)CCNeighbors[i][j] << " " << coildist << endl;
           }
      }

	}
  }

	outfile.close();
}



