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
void lpGenerator::writeLP(string alphaTSname, string betaTSname, string coilTSname, ofstream* output){

  //createMapping(lengthOfWord);
  cout << "I'm at the beginning of LP writing" << endl;

  writeObjectiveFunction(alphaTSname, betaTSname, coilTSname, output, TRSSize);
  sleep(5 * 60);

  writeSigmaFunctions(output);
  sleep(5 * 60);
  writeErrorFunctionsFromFiles(output);
  //deleteMapping(lengthOfWord);
}

void lpGenerator::buildSaveTrees(char* alpha, char* beta, char* alphaSave, char* betaSave,int rand,bool weights) {
  const char* AlphaInputFile = (const char* ) alpha;
  const char* BetaInputFile = (const char* ) beta;
  string r = "r";
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  if (weights){
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
      thread alphatree(&lpGenerator::buildalphatreeweighted, this, alphaFile, 1, false, rand);
      thread betatree(&lpGenerator::buildbetatreeweighted, this, betaFile, 1, false, rand);
      alphatree.join();
      betatree.join();
  }else{
      thread alphatree(&lpGenerator::buildalphatree, this, alphaFile, 1, false, rand);
      thread betatree(&lpGenerator::buildbetatree, this, betaFile, 1, false, rand);
      alphatree.join();
      betatree.join();
  }
  nnfinder1.saveDispersionTree(alphaSave, "alpha");
  nnfinder4.saveDispersionTree(betaSave, "beta");
}
void lpGenerator::buildSaveTreesCoil(char* alpha, char* beta,char* coil, char* alphaSave, char* betaSave,char* coilSave,int rand,bool weights) {
  const char* AlphaInputFile = (const char* ) alpha;
  const char* BetaInputFile = (const char* ) beta;
  const char* CoilInputFile = (const char* ) coil;
  string r = "r";
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  FILE* coilFile = fopen(CoilInputFile, r.c_str());
  if (weights){
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
      thread alphatree(&lpGenerator::buildalphatreeweighted, this, alphaFile, 1, false, rand);
      thread betatree(&lpGenerator::buildbetatreeweighted, this, betaFile, 1, false, rand);
      thread coiltree(&lpGenerator::buildcoiltreeweighted, this, coilFile, 1, false, rand);
      alphatree.join();
      betatree.join();
      coiltree.join();
  }else{
      thread alphatree(&lpGenerator::buildalphatree, this, alphaFile, 1, false, rand);
      thread betatree(&lpGenerator::buildbetatree, this, betaFile, 1, false, rand);
      thread coiltree(&lpGenerator::buildcoiltree, this, coilFile, 1, false, rand);
      alphatree.join();
      betatree.join();
      coiltree.join();
  }
  nnfinder1.saveDispersionTree(alphaSave, "alpha");
  nnfinder4.saveDispersionTree(betaSave, "beta");
  nnfinder7.saveDispersionTree(coilSave, "coil");
}
void lpGenerator::buildSStrees(char* alpha, char* beta,char* coil, char* alphaSave, char* betaSave,char* coilSave,int rand,bool weights) {
  const char* AlphaInputFile = (const char* ) alpha;
  const char* BetaInputFile = (const char* ) beta;
  const char* CoilInputFile = (const char* ) coil;
  string r = "r";
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  FILE* coilFile = fopen(CoilInputFile, r.c_str());
  if (weights){
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
      thread alphatree(&lpGenerator::buildcombinedalphatree, this, alphaFile, 1, false, rand);
      thread betatree(&lpGenerator::buildcombinedbetatree, this, betaFile, 1, false, rand);
      thread coiltree(&lpGenerator::buildcombinedcoiltree, this, coilFile, 1, false, rand);
      alphatree.join();
      betatree.join();
      coiltree.join();
  }
  else {
      thread alphatree(&lpGenerator::buildcombinedalphatree, this, alphaFile, 1, false, rand);
      thread betatree(&lpGenerator::buildcombinedbetatree, this, betaFile, 1, false, rand);
      thread coiltree(&lpGenerator::buildcombinedcoiltree, this, coilFile, 1, false, rand);
      alphatree.join();
      betatree.join();
      coiltree.join();
  }
  nnfinder1.saveDispersionTree(alphaSave, "alpha");
  nnfinder4.saveDispersionTree(betaSave, "beta");
  nnfinder7.saveDispersionTree(coilSave, "coil");
}
void lpGenerator::printDFs(){
    //nnfinder1.printDFs();
}
void lpGenerator::buildSaveTreesDefault(char* alpha, char* beta, char* alphaSave, char* betaSave,int rand,bool weights) {
  const char* AlphaInputFile = (const char* ) alpha;
  const char* BetaInputFile = (const char* ) beta;
  string r = "r";
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  setEpsilonRho(1,0.9,0.3);
  if (weights){
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
      thread alphatree(&lpGenerator::buildalphatreeweighted, this, alphaFile, 1, true, rand);
      thread betatree(&lpGenerator::buildbetatreeweighted, this, betaFile, 1, true, rand);
      alphatree.join();
      betatree.join();
  }else{
      thread alphatree(&lpGenerator::buildalphatree, this, alphaFile, 1, true, rand);
      thread betatree(&lpGenerator::buildbetatree, this, betaFile, 1, true, rand);
      alphatree.join();
      betatree.join();
  }
  nnfinder1.saveDispersionTree(alphaSave, "alpha");
  nnfinder4.saveDispersionTree(betaSave, "beta");
}
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
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
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

void lpGenerator::buildmap(string TrainingProteins){
    ifstream trainingfile(TrainingProteins.c_str());
    int i = 0;
    string line;
    unsigned int wordlen = 23;
    cout << "starting building map" << endl;
    while( !trainingfile.eof()){
        getline(trainingfile,line);
        if (trainingfile){
            ++i;
            if (i % 5 == 2 && line.length() >= wordlen){
                //process stuff
                for (unsigned int j = 0;j <= line.length() - wordlen;j++){
                    string word = line.substr(j,wordlen);
                    if(this->countmap.find(word) == this->countmap.end()){
                        this->countmap[word] = 1;
                    }
                    else{
                        this->countmap[word] = this->countmap[word] + 1;
                    }
                }
            }
        }
        else break;
    }
}



void lpGenerator::findNeighborsOneClass(char* TreeIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand, float rho, bool counts, bool load) {
  const char* TreeFilename = (const char* ) TreeIF;
  const char* QueryFilename = (const char* ) QueryIF;
  float epsilon = 1;
  setEpsilonRho(epsilon,rho,0.3);
  //cout << "rhos in use are: aarho = " << aarho << " ssrho = " << ssrho << endl;
  string r = "r";
  FILE* TreeFile = fopen(TreeFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());

  if (load == true){
      cout << "loading trees" << endl;
      if(firstIteration == true){
          if (type == 1)
              nnfinder1.loadDispersionTree(TreeFile, "trees/alpha" + to_string(rand) + "4.tree", "defaultaw",rand);
          else if (type == 2)
              nnfinder4.loadDispersionTree(TreeFile, "trees/beta" + to_string(rand) + "4.tree","defaultbw",rand);
          else
              nnfinder7.loadDispersionTree(TreeFile, "trees/coil" + to_string(rand) + "4.tree","defaultcw",rand);

      }
      else{
          if (type == 1)
              nnfinder1.loadDispersionTree(TreeFile, "trees/alpha" + to_string(rand) + ".tree", "alpha",rand);
          else if (type == 2)
              nnfinder4.loadDispersionTree(TreeFile, "trees/beta" + to_string(rand) + ".tree","beta",rand);
          else
              nnfinder7.loadDispersionTree(TreeFile, "trees/coil" + to_string(rand) + ".tree","coil",rand);

      }

  }
  else{

      cout << "building trees" << endl;
      if (type == 1)
          this->buildalphatree( TreeFile, 1, firstIteration, rand);
      else
          this->buildbetatree(TreeFile, 1, firstIteration, rand);
  }


  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  MetricSpace* MSQuery = CreateMetricSpace(inputchar, QueryFile);
  lpGenerator::TRSSize = MSQuery->numpoints;
  cout << "built the trees" << " and numpoints = " << MSQuery->numpoints << endl;

  //thread t1(&lpGenerator::threadalphafunc, this, k + 1, MSAlpha, &this->AANeighbors, 1);
  //thread t2(&lpGenerator::threadalphafunc, this,l, MSBeta2, &this->ABNeighbors, 2);
  //thread t3(&lpGenerator::threadbetafunc, this,k + 1, MSBeta, &this->BBNeighbors, 1);
  //thread t4(&lpGenerator::threadbetafunc, this, l, MSAlpha2, &this->BANeighbors, 2);
  //thread t5(&lpGenerator::threadalphafunc, this,l, MSCoil, &this->ACNeighbors, 3);
  //thread t6(&lpGenerator::threadbetafunc, this,l, MSCoil2, &this->BCNeighbors, 3);

  //t1.join();
  //t2.join();
  //t3.join();
  //t4.join();
  //t5.join();
  //t6.join();
  if (type == 1){
      nnfinder1.queryAlpha(numneighbors + 1, MSQuery, &this->AANeighbors, 1);
      
  }
  else if (type == 2){
      nnfinder4.queryBeta(numneighbors + 1, MSQuery, &this->AANeighbors, 1);
  }
  else if (type == 3){
      nnfinder7.queryCoil(numneighbors + 1, MSQuery, &this->AANeighbors, 1);
  }
  cout << "all threads joined" << endl;
  ofstream outfile(outputfile.c_str());
  ifstream infile(QueryFilename);
  string inword;
  infile >> inword;
  infile >> inword;
  infile >> inword;
  for (int i = 0;i < TRSSize;i++) {
      cout << i << endl;
    infile >> inword;
    outfile << inword << " ";
    int numprinted = 0;
    int j = 0;
    int curcount = 0;
    if(counts){
        cout << "j " << j <<endl;
        string q((char*)this->AANeighbors[i][j]);
        if(this->countmap.find(q) != this->countmap.end()){
            curcount = this->countmap[q];
        }
        else{
            curcount = 1;
        }

        if(strcmp(inword.c_str(),(char*) this->AANeighbors[i][j]) == 0)
            curcount--;
    }
    while(numprinted < numneighbors){
        if (counts){
            if (curcount > 0){
                outfile << (char*) this->AANeighbors[i][j];
                numprinted++;
                if (numprinted != numneighbors)
                    outfile << " ";
                curcount = 0;
            }
            else{
                j++;
                string q((char*)this->AANeighbors[i][j]);
                if(this->countmap.find(q) != this->countmap.end()){
                    curcount = this->countmap[q];
                }
                else{
                    curcount = 1;
                }
                if(strcmp(inword.c_str(),(char*) this->AANeighbors[i][j]) == 0)
                    curcount--;
            }

        }
        else{
            if(strcmp(inword.c_str(),(char*) this->AANeighbors[i][j]) != 0){
              outfile << (char*) this->AANeighbors[i][j];
              numprinted++;
              if (numprinted != numneighbors)
                  outfile << " ";
            }
            j++;
        }
    }
    outfile << endl;
  }
  outfile.close();
  infile.close();

  DestroyMetricSpace(MSQuery);
  fclose(TreeFile);
  fclose(QueryFile);
  delete[] inputchar;
}

void lpGenerator::findNeighborsOneClassCombined(char* TreeIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand, float epsilon, float aarho, float ssrho) {
  const char* TreeFilename = (const char* ) TreeIF;
  const char* QueryFilename = (const char* ) QueryIF;
  setEpsilonRho(epsilon,aarho,ssrho);
  string r = "r";
  FILE* TreeFile = fopen(TreeFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());
  cout << "printing the wordlength:" << lengthOfWord <<endl;

  if (numneighbors < 10){
      if (type == 1){

          this->buildcombinedalphatree( TreeFile, 1, firstIteration,rand);
      }
      else{
          this->buildcombinedbetatree( TreeFile, 1, firstIteration,rand);
      }
  }
  else{
      if (type == 1){
          this->buildalphatree( TreeFile,1,firstIteration,rand);
      }
      else{
          this->buildbetatree( TreeFile,1,firstIteration,rand);
      }
  }


  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  MetricSpace* MSQuery = CreateMetricSpace(inputchar, QueryFile);
  //cout
  lpGenerator::TRSSize = MSQuery->numpoints;
  cout << "built the trees" << endl;

  //thread t1(&lpGenerator::threadalphafunc, this, k + 1, MSAlpha, &this->AANeighbors, 1);
  //thread t2(&lpGenerator::threadalphafunc, this,l, MSBeta2, &this->ABNeighbors, 2);
  //thread t3(&lpGenerator::threadbetafunc, this,k + 1, MSBeta, &this->BBNeighbors, 1);
  //thread t4(&lpGenerator::threadbetafunc, this, l, MSAlpha2, &this->BANeighbors, 2);
  //thread t5(&lpGenerator::threadalphafunc, this,l, MSCoil, &this->ACNeighbors, 3);
  //thread t6(&lpGenerator::threadbetafunc, this,l, MSCoil2, &this->BCNeighbors, 3);

  //t1.join();
  //t2.join();
  //t3.join();
  //t4.join();
  //t5.join();
  //t6.join();
  if (type == 1){
      nnfinder1.queryAlpha(numneighbors+1, MSQuery, &this->AANeighbors, 1);
      
  }
  else{
      nnfinder4.queryBeta(numneighbors+1, MSQuery, &this->AANeighbors, 1);
  }
  cout << "all threads joined" << endl;
  ofstream outfile(outputfile.c_str());
  ifstream infile(QueryFilename);
  string inword;
  infile >> inword;
  infile >> inword;
  infile >> inword;
  char* outword = new char[lengthOfWord + 1];
  char* queryword = new char[lengthOfWord + 1];
  for (int i = 0;i < TRSSize;i++) {
    infile >> inword;
    strncpy(queryword,inword.c_str(),lengthOfWord);
    outfile << queryword << " ";

    int numwritten = 0;
    int j = 0;
    while (numwritten < numneighbors) {
        char* combinedword = (char*) this->AANeighbors[i][j];
        strncpy(outword,combinedword,lengthOfWord);
        outword[lengthOfWord] = '\0';
        if (strcmp(outword,queryword) != 0 && numwritten < numneighbors){
          outfile << outword;
          numwritten++;
        
      if (j != numneighbors - 1)
          outfile << " ";
        }
        j++;
    }
    //infile >> inword;
    outfile << endl;
  }
  outfile.close();
  infile.close();

  delete[] outword;
  DestroyMetricSpace(MSQuery);
  fclose(TreeFile);
  fclose(QueryFile);
  delete[] inputchar;
}
void lpGenerator::findNeighborsCoverset(char* TreeIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand) {
  const char* TreeFilename = (const char* ) TreeIF;
  const char* QueryFilename = (const char* ) QueryIF;
  string r = "r";
  FILE* TreeFile = fopen(TreeFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());
  setEpsilonRho(1,0.9,0.3);

  if (type == 1){

      this->buildalphatree( TreeFile, 1, firstIteration,rand);
  }
  else{
      this->buildbetatree( TreeFile, 1, firstIteration,rand);
  }


  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  MetricSpace* MSQuery = CreateMetricSpace(inputchar, QueryFile);
  lpGenerator::TRSSize = MSQuery->numpoints;
  cout << "built the trees" << endl;

  //thread t1(&lpGenerator::threadalphafunc, this, k + 1, MSAlpha, &this->AANeighbors, 1);
  //thread t2(&lpGenerator::threadalphafunc, this,l, MSBeta2, &this->ABNeighbors, 2);
  //thread t3(&lpGenerator::threadbetafunc, this,k + 1, MSBeta, &this->BBNeighbors, 1);
  //thread t4(&lpGenerator::threadbetafunc, this, l, MSAlpha2, &this->BANeighbors, 2);
  //thread t5(&lpGenerator::threadalphafunc, this,l, MSCoil, &this->ACNeighbors, 3);
  //thread t6(&lpGenerator::threadbetafunc, this,l, MSCoil2, &this->BCNeighbors, 3);

  //t1.join();
  //t2.join();
  //t3.join();
  //t4.join();
  //t5.join();
  //t6.join();

  if (type == 1){
      nnfinder1.queryAlpha(numneighbors, MSQuery, &this->AANeighbors, 1);
      
  }
  else{
      nnfinder4.queryBeta(numneighbors, MSQuery, &this->AANeighbors, 1);
  }
  cout << "all threads joined" << endl;
  ofstream outfile(outputfile.c_str());
  ifstream infile(QueryFilename);
  string inword;
  infile >> inword;
  infile >> inword;
  infile >> inword;
  for (int i = 0;i < TRSSize;i++) {
    infile >> inword;
    outfile << inword << " ";
    if(firstIteration == true){
        for (int j = 0;j < numneighbors; j++) {
          outfile << (char*) this->AANeighbors[i][j] << " " << WordDist(MSQuery->points[i],this->AANeighbors[i][j]);
          if (j != numneighbors - 1)
              outfile << " ";
        }
    }
    else{
        if(type == 1){
            for (int j = 0;j < numneighbors; j++) {
              outfile << (char*) this->AANeighbors[i][j] << " " << wordDist2Alpha(MSQuery->points[i],this->AANeighbors[i][j]);
              if (j != numneighbors - 1)
                  outfile << " ";
            }
        }
        else{
            for (int j = 0;j < numneighbors; j++) {
              outfile << (char*) this->AANeighbors[i][j] << " " << wordDist2Beta(MSQuery->points[i],this->AANeighbors[i][j]);
              if (j != numneighbors - 1)
                  outfile << " ";
            }
        }
    }

    outfile << endl;
  }
  outfile.close();
  infile.close();

  DestroyMetricSpace(MSQuery);
  fclose(TreeFile);
  fclose(QueryFile);
  delete[] inputchar;
}

void lpGenerator::findDistances2(char* QF, char* NF, string outputfile, int rand){

  const char* NeighborsFilename = (const char* ) NF;
  const char* QueryFilename = (const char* ) QF;
  string r = "r";
  FILE* NeighborsFile = fopen(NeighborsFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());
  string alphadfstring = "dfs/alphadis" + to_string(rand) + ".h";
  string betadfstring = "dfs/betadis" + to_string(rand) + ".h";
  char* adf = new char[alphadfstring.size() + 1];
  copy(alphadfstring.begin(), alphadfstring.end(), adf);
  adf[alphadfstring.size()] = '\0';
  char* bdf = new char[betadfstring.size() + 1];
  copy(betadfstring.begin(), betadfstring.end(), bdf);
  bdf[betadfstring.size()] = '\0';
  //alphadf = adf;
  //betadf = bdf;
  cout << endl << "setting distance functions with rand = " << rand << endl;
  setdfs(adf,bdf);

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  MetricSpace* MSQuery = CreateMetricSpace(inputchar, QueryFile);
  MetricSpace* Neighbors = CreateMetricSpace(inputchar, NeighborsFile);
  ofstream outfile(outputfile.c_str());
  for (unsigned int i = 0;i < MSQuery->numpoints;i++) {
    outfile << wordDist2Alpha(MSQuery->points[i],Neighbors->points[i]) << " " << wordDist2Beta(MSQuery->points[i],Neighbors->points[i]) << " " ;
    outfile << endl;
    
    }

    outfile.close();
    DestroyMetricSpace(MSQuery);
    DestroyMetricSpace(Neighbors);
    fclose(QueryFile);
    fclose(NeighborsFile);
    delete[] inputchar;



    



}

void lpGenerator::findCombinedDistances2(char* QF, char* NF, string outputfile, int rand,float epsilon, float aarho, float ssrho){

    setEpsilonRho(epsilon,aarho, ssrho);
  const char* NeighborsFilename = (const char* ) NF;
  const char* QueryFilename = (const char* ) QF;
  string r = "r";
  FILE* NeighborsFile = fopen(NeighborsFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

    string alphadfstring = "dfs/alphadis" + to_string(rand) + ".h";
    string betadfstring = "dfs/betadis" + to_string(rand) + ".h";
    char* adf = new char[alphadfstring.size() + 1];
    copy(alphadfstring.begin(), alphadfstring.end(), adf);
    adf[alphadfstring.size()] = '\0';
    char* bdf = new char[betadfstring.size() + 1];
    copy(betadfstring.begin(), betadfstring.end(), bdf);
    bdf[betadfstring.size()] = '\0';
    //alphadf = adf;
    //betadf = bdf;
    setdfs(adf,bdf);
     
  MetricSpace* MSQuery = CreateMetricSpace(inputchar, QueryFile);
  MetricSpace* Neighbors = CreateMetricSpace(inputchar, NeighborsFile);
  ofstream outfile(outputfile.c_str());
  for (unsigned int i = 0;i < MSQuery->numpoints;i++) {
    outfile << CombinedAlphaDistance(MSQuery->points[i],Neighbors->points[i]) << " " << CombinedBetaDistance(MSQuery->points[i],Neighbors->points[i]) << " " << CombinedWordDist(MSQuery->points[i],Neighbors->points[i]) << " ";
    outfile << endl;
    
    }

    outfile.close();
    DestroyMetricSpace(MSQuery);
    DestroyMetricSpace(Neighbors);
    fclose(QueryFile);
    fclose(NeighborsFile);
    delete[] inputchar;



    



}
void lpGenerator::findDistances(string QueryFilename,string outputfile)
{
    
  ifstream QueryFile(QueryFilename.c_str());
  string** inputwords = NULL;
  string line;
  int wordListSize = 0;
  string firstline;
  getline(QueryFile,firstline);
  cout << firstline << endl;
  //firstline.erase(firstline.find("\n"));
  int sizeofline = 0;
  while(getline(QueryFile,line)){
      sizeofline = 1;
      size_t next = line.find(" ");
      size_t last = 0;
      while(next != string::npos){
          sizeofline++;
          last = next;
          cout << last << endl;
          next = line.find(" ",last + 1);
          cout << next << endl;
      }

      inputwords[wordListSize] = new string[sizeofline];
      cout << line << endl;
      string firstword = line.substr(0,line.find(" "));
      inputwords[wordListSize][0] = firstword;
      line = line.erase(0,line.find(" ") + 1);
      int rowlen = 1;
      while(line.find(" ") != string::npos){
          string nextword = line.substr(0,line.find(" "));
          inputwords[wordListSize][rowlen] = nextword;
          rowlen++;
      }
      wordListSize++;
  }

  ofstream outfile(outputfile.c_str());
  for (int i = 0;i < wordListSize;i++) {
      for (int j = 1; j < sizeofline;j++){
        outfile << wordDist2Alpha((Pointer)inputwords[i][0].c_str(),(Pointer)inputwords[i][j].c_str()) << " " << wordDist2Beta((Pointer)inputwords[i][0].c_str(),(Pointer)inputwords[i][j].c_str()) << " " ;
        outfile << "one";
    
    }
    outfile << endl;

  }
  for (int i = 0; i < wordListSize;i++){
      delete[] inputwords[i];
  }
  QueryFile.close();
  outfile.close();
}

void lpGenerator::findNeighborsOneClassDistance(char* TreeIF,char* BetaIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand) {
  const char* TreeFilename = (const char* ) TreeIF;
  const char* betaTreeFilename = (const char*) BetaIF;
  const char* QueryFilename = (const char* ) QueryIF;
  string r = "r";
  FILE* TreeFile = fopen(TreeFilename, r.c_str());
  FILE* BetaFile = fopen(betaTreeFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());

  this->buildalphatree( TreeFile, 1, firstIteration, rand);
  this->buildbetatree(BetaFile, 1, firstIteration, rand);

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  MetricSpace* MSQuery = CreateMetricSpace(inputchar, QueryFile);
  lpGenerator::TRSSize = MSQuery->numpoints;
  cout << "built the trees" << endl;

  //thread t1(&lpGenerator::threadalphafunc, this, k + 1, MSAlpha, &this->AANeighbors, 1);
  //thread t2(&lpGenerator::threadalphafunc, this,l, MSBeta2, &this->ABNeighbors, 2);
  //thread t3(&lpGenerator::threadbetafunc, this,k + 1, MSBeta, &this->BBNeighbors, 1);
  //thread t4(&lpGenerator::threadbetafunc, this, l, MSAlpha2, &this->BANeighbors, 2);
  //thread t5(&lpGenerator::threadalphafunc, this,l, MSCoil, &this->ACNeighbors, 3);
  //thread t6(&lpGenerator::threadbetafunc, this,l, MSCoil2, &this->BCNeighbors, 3);

  //t1.join();
  //t2.join();
  //t3.join();
  //t4.join();
  //t5.join();
  //t6.join();
  nnfinder1.queryAlpha(numneighbors, MSQuery, &this->AANeighbors, 1);
  nnfinder4.queryBeta(numneighbors, MSQuery, &this->BBNeighbors, 1);
  cout << "all threads joined" << endl;
  ofstream outfile(outputfile.c_str());
  for (int i = 0;i < TRSSize;i++) {
	  if ((type == 2 && wordDist2Beta(MSQuery->points[i], this->BBNeighbors[i][1]) > wordDist2Alpha(MSQuery->points[i], this->AANeighbors[i][0])) || (type == 1 && wordDist2Beta(MSQuery->points[i], this->BBNeighbors[i][0]) > wordDist2Alpha(MSQuery->points[i], this->AANeighbors[i][1])) || (type == 3 && wordDist2Beta(MSQuery->points[i], this->BBNeighbors[i][0]) < wordDist2Alpha(MSQuery->points[i], this->AANeighbors[i][0]))  ){
  	if (type != 1)
	  outfile << (char*) MSQuery->points[i] << " ";
    for (int j = 0;j < numneighbors; j++) {
	    if (type == 1 || j == 1)
      		outfile << (char*) this->AANeighbors[i][j] << " ";
    }
    if (type == 1)
    	outfile <<  wordDist2Alpha(this->AANeighbors[i][0], this->AANeighbors[i][1]);
    else
	    outfile << wordDist2Alpha(MSQuery->points[i], this->AANeighbors[i][0]);
    outfile << endl;
    }
  }
  outfile.close();

  DestroyMetricSpace(MSQuery);
  fclose(TreeFile);
  fclose(QueryFile);
  delete[] inputchar;
}

void lpGenerator::findNeighborsOneClassTwoDistance(char* AlphaIF, char* BetaIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand, bool load, bool weights, bool training) {
  const char* AlphaFilename = (const char* ) AlphaIF;
  const char* BetaFilename = (const char* ) BetaIF;
  const char* QueryFilename = (const char* ) QueryIF;
  string r = "r";
  FILE* AlphaFile = fopen(AlphaFilename, r.c_str());
  FILE* BetaFile = fopen(BetaFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());

  if (firstIteration == true){
      float epsilon = 1;
      float aarho = 0.9;
      float ssrho = 0.3;
      setEpsilonRho(epsilon,aarho, ssrho);
  }
  if(weights == true){
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
  }
  if (load == true){
      cout << "loading trees" << endl;
      if(firstIteration == true){
          if (weights){
              nnfinder1.loadDispersionTree(AlphaFile, "trees/alpha"+ to_string(rand) + "4.tree", "defaultaw",rand);
              nnfinder4.loadDispersionTree(BetaFile, "trees/beta" + to_string(rand) + "4.tree","defaultbw",rand);
          }
          else{
              nnfinder1.loadDispersionTree(AlphaFile, "trees/alpha04.tree", "defaulta",rand);
              nnfinder4.loadDispersionTree(BetaFile, "trees/beta04.tree","defaultb",rand);
          }
      }
      else if(weights == true){
          nnfinder1.loadDispersionTree(AlphaFile, "trees/alpha" + to_string(rand) + "4.tree","alphaweighted",rand);
          nnfinder4.loadDispersionTree(BetaFile, "trees/beta" + to_string(rand) + "4.tree","betaweighted",rand);
      }
      else{
          nnfinder1.loadDispersionTree(AlphaFile, "trees/alpha" + to_string(rand) + "4.tree", "alpha",rand);
          nnfinder4.loadDispersionTree(BetaFile, "trees/beta" + to_string(rand) + "4.tree","beta",rand);

      }

  }
  else{

      if(weights == true){
          cout << "building weights trees" << endl;
          this->buildalphatreeweighted(AlphaFile,1,firstIteration,rand);
          this->buildbetatreeweighted(BetaFile,1,firstIteration,rand);

      }
      else{
          cout << "building trees" << endl;
          this->buildalphatree( AlphaFile, 1, firstIteration, rand);
          this->buildbetatree(BetaFile, 1, firstIteration, rand);
      }
  }
  Distance alphadistfunc;
  Distance betadistfunc;
  if(firstIteration == true){
      alphadistfunc = WordDist;
      betadistfunc = WordDist;
  }
  else if(weights == true){
      alphadistfunc = WeightedAlphaDistance;
      betadistfunc = WeightedBetaDistance;
  }
  else{
      alphadistfunc = wordDist2Alpha;
      betadistfunc = wordDist2Beta;
  }
  cout << "I'm finding the neighbors and their two distances" << endl;

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  MetricSpace* MSQuery = CreateMetricSpace(inputchar, QueryFile);
  lpGenerator::TRSSize = MSQuery->numpoints;
  cout << "built the trees" << endl;

  nnfinder1.queryAlpha(numneighbors, MSQuery, &this->AANeighbors, 1);
  nnfinder4.queryBeta(numneighbors, MSQuery, &this->BBNeighbors, 1);
  cout << "all threads joined" << endl;
  string a = outputfile + "a";
  ofstream outfilea(a.c_str());
  for (int i = 0;i < TRSSize;i++) {
      for (int j = 0; j < numneighbors - 1; j++){
          if (type == 1){
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->AANeighbors[i][j + 1] << " " << alphadistfunc(MSQuery->points[i], this->AANeighbors[i][j + 1]) << " " ;
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->BBNeighbors[i][j] << " " << betadistfunc(MSQuery->points[i], this->BBNeighbors[i][j]) << endl;
            }
          else if (type == 2) {
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->AANeighbors[i][j] << " " << alphadistfunc(MSQuery->points[i], this->AANeighbors[i][j]) << " " ;
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->BBNeighbors[i][j+1] << " " << betadistfunc(MSQuery->points[i], this->BBNeighbors[i][j+1]) << endl ;
          }
          else if (type == 3) {
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->AANeighbors[i][j] << " " << alphadistfunc(MSQuery->points[i], this->AANeighbors[i][j]) << " " ;
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->BBNeighbors[i][j] << " " << betadistfunc(MSQuery->points[i], this->BBNeighbors[i][j]) << endl;
          }
      }
  }
  outfilea.close();

  DestroyMetricSpace(MSQuery);
  fclose(AlphaFile);
  fclose(BetaFile);
  fclose(QueryFile);
  delete[] inputchar;
}
void lpGenerator::findNeighborsCoilClassTwoDistance(char* CoilIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand, bool load, bool weights, bool training) {
  const char* CoilFilename = (const char* ) CoilIF;
  const char* QueryFilename = (const char* ) QueryIF;
  string r = "r";
  FILE* CoilFile = fopen(CoilFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());

  if (firstIteration == true){
      float epsilon = 1;
      float aarho = 0.9;
      float ssrho = 0.3;
      setEpsilonRho(epsilon,aarho, ssrho);
  }
  if(weights == true){
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
  }
  if (load == true){
      cout << "loading trees" << endl;
      if(firstIteration == true){
          if (weights){
              nnfinder7.loadDispersionTree(CoilFile, "trees/coil"+ to_string(rand) + "4.tree", "defaultcw",rand);
          }
          else{
              nnfinder7.loadDispersionTree(CoilFile, "trees/coil04.tree", "defaultc",rand);
          }
      }
      else if(weights == true){
          nnfinder7.loadDispersionTree(CoilFile, "trees/coil" + to_string(rand) + "4.tree","coilweighted",rand);
      }
      else{
          nnfinder7.loadDispersionTree(CoilFile, "trees/coil" + to_string(rand) + "4.tree", "coil",rand);

      }

  }
  else{

      if(weights == true){
          cout << "building weights trees" << endl;
          this->buildcoiltreeweighted(CoilFile,1,firstIteration,rand);

      }
      else{
          cout << "building trees" << endl;
          this->buildcoiltree( CoilFile, 1, firstIteration, rand);
      }
  }
  Distance coildistfunc;
  if(firstIteration == true){
      coildistfunc = WordDist;
  }
  else if(weights == true){
      coildistfunc = WeightedCoilDistance;
  }
  else{
      coildistfunc = wordDist2Coil;
  }
  cout << "I'm finding the neighbors and their two distances" << endl;

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  MetricSpace* MSQuery = CreateMetricSpace(inputchar, QueryFile);
  lpGenerator::TRSSize = MSQuery->numpoints;
  cout << "built the trees" << endl;

  nnfinder7.queryCoil(numneighbors, MSQuery, &this->CCNeighbors, 1);
  cout << "all threads joined" << endl;
  string a = outputfile + "a";
  ofstream outfilea(a.c_str());
  for (int i = 0;i < TRSSize;i++) {
      for (int j = 0; j < numneighbors - 1; j++){
          if (type == 1){
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->CCNeighbors[i][j] << " " << coildistfunc(MSQuery->points[i], this->CCNeighbors[i][j + 1]) << endl ;
            }
          else if (type == 2) {
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->CCNeighbors[i][j] << " " << coildistfunc(MSQuery->points[i], this->CCNeighbors[i][j]) << endl ;
          }
          else if (type == 3) {
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->CCNeighbors[i][j + 1] << " " << coildistfunc(MSQuery->points[i], this->CCNeighbors[i][j]) << endl ;
          }
      }
  }
  outfilea.close();

  DestroyMetricSpace(MSQuery);
  fclose(CoilFile);
  fclose(QueryFile);
  delete[] inputchar;
}

 
void lpGenerator::findNeighborsOneClassTwoDistanceAve(char* AlphaIF, char* BetaIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand, bool counts) {
  const char* AlphaFilename = (const char* ) AlphaIF;
  const char* BetaFilename = (const char* ) BetaIF;
  const char* QueryFilename = (const char* ) QueryIF;
  string r = "r";
  FILE* AlphaFile = fopen(AlphaFilename, r.c_str());
  FILE* BetaFile = fopen(BetaFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());
  if (firstIteration == true){
      float epsilon = 1;
      float aarho = 0.9;
      float ssrho = 0.3;
      setEpsilonRho(epsilon,aarho, ssrho);
  }

  cout << "building trees" << endl;
  this->buildalphatree( AlphaFile, 1, firstIteration, rand);
  this->buildbetatree(BetaFile, 1, firstIteration, rand);
  cout << "I'm finding the neighbors and their two distances" << endl;

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  MetricSpace* MSQuery = CreateMetricSpace(inputchar, QueryFile);
  lpGenerator::TRSSize = MSQuery->numpoints;
  cout << "built the trees" << endl;

  nnfinder1.queryAlpha(numneighbors, MSQuery, &this->AANeighbors, 1);
  nnfinder4.queryBeta(numneighbors, MSQuery, &this->BBNeighbors, 1);
  cout << "all threads joined" << endl;
  string a = outputfile + "a";
  ofstream outfilea(a.c_str());
  for (int i = 0;i < TRSSize;i++) {
      float alphadist = 0;
      float betadist = 0;
      string inword( (char*) MSQuery->points[i]);
      if(counts){
          int anumsummed = 0;
          int aj =0;
          int acurcount = 0;
          int bnumsummed = 0;
          int bj =0;
          int bcurcount = 0;
        string q((char*)this->AANeighbors[i][aj]);
        if(this->countmap.find(q) != this->countmap.end()){
            acurcount = this->countmap[q];
        }
        else{
            acurcount = 1;
            cout << q << " was not in the map but was in the touchstone" << endl;
        }

        if(strcmp(inword.c_str(),(char*) this->AANeighbors[i][aj]) == 0)
            acurcount--;
        string p((char*)this->BBNeighbors[i][bj]);
        if(this->countmap.find(p) != this->countmap.end()){
            bcurcount = this->countmap[p];
        }
        else{
            bcurcount = 1;
            cout << p << " was not in the map but was in the touchstone" << endl;
        }

        if(strcmp(inword.c_str(),(char*) this->BBNeighbors[i][bj]) == 0)
            bcurcount--;
          while(anumsummed < numneighbors || bnumsummed < numneighbors){
              if (acurcount > 0 && anumsummed < numneighbors){
                  if (firstIteration)
                      alphadist += WordDist(MSQuery->points[i],this->AANeighbors[i][aj]);
                  else
                      alphadist += wordDist2Alpha(MSQuery->points[i],this->AANeighbors[i][aj]);
                  acurcount--;
                  anumsummed++;
              }
              else{
                  aj++;
                  string u((char*)this->AANeighbors[i][aj]);
                  if(this->countmap.find(u) != this->countmap.end()){
                      acurcount = this->countmap[u];
                  }
                  else{
                      acurcount = 1;
                      cout << u << " was not in the map but was in the touchstone" << endl;
                  }

                  if(strcmp(inword.c_str(),(char*) this->AANeighbors[i][aj]) == 0)
                      acurcount--;

              }
              if (bcurcount > 0 && bnumsummed < numneighbors){
                  if (firstIteration)
                      betadist += WordDist(MSQuery->points[i],this->BBNeighbors[i][bj]);
                  else
                      betadist += wordDist2Beta(MSQuery->points[i],this->BBNeighbors[i][bj]);
                  bcurcount--;
                  bnumsummed++;
              }
              else{
                  bj++;
                  string t((char*)this->BBNeighbors[i][bj]);
                  if(this->countmap.find(t) != this->countmap.end()){
                      bcurcount = this->countmap[t];
                  }
                  else{
                      bcurcount = 1;
                      cout << t << " was not in the map but was in the touchstone" << endl;
                  }
  
                  if(strcmp(inword.c_str(),(char*) this->BBNeighbors[i][bj]) == 0)
                      bcurcount--;
              }
          }

      }else{
          for (int j = 0; j < numneighbors - 1; j++){
              if (type == 1){
                  if(firstIteration == false){
                      alphadist += wordDist2Alpha(MSQuery->points[i], this->AANeighbors[i][j+1]);
                      betadist += wordDist2Beta(MSQuery->points[i], this->BBNeighbors[i][j]);
                  }
                  else{
                      alphadist += WordDist(MSQuery->points[i], this->AANeighbors[i][j+1]);
                      betadist += WordDist(MSQuery->points[i], this->BBNeighbors[i][j]);
                  }

                }
              else if (type == 2) {
                  if(firstIteration == false){
                      alphadist += wordDist2Alpha(MSQuery->points[i], this->AANeighbors[i][j]);
                      betadist += wordDist2Beta(MSQuery->points[i], this->BBNeighbors[i][j+1]);
                  }
                  else{
                      alphadist += WordDist(MSQuery->points[i], this->AANeighbors[i][j]);
                      betadist += WordDist(MSQuery->points[i], this->BBNeighbors[i][j+1]);
                  }

              }
              else if (type == 3) {
                  if(firstIteration == false){
                      alphadist += wordDist2Alpha(MSQuery->points[i], this->AANeighbors[i][j]);
                      betadist += wordDist2Beta(MSQuery->points[i], this->BBNeighbors[i][j]);
                  }
                  else{
                      alphadist += WordDist(MSQuery->points[i], this->AANeighbors[i][j]);
                      betadist += WordDist(MSQuery->points[i], this->BBNeighbors[i][j]);
                  }
              }
          }
      }
      outfilea << (char*) MSQuery->points[i] << " ";
      if (counts)
          outfilea <<(char*) this->AANeighbors[i][1] << " " << alphadist / (1.0 * (numneighbors)) << " " ;
      else
          outfilea <<(char*) this->AANeighbors[i][1] << " " << alphadist / (1.0 * (numneighbors-1)) << " " ;
      outfilea << (char*) MSQuery->points[i] << " ";
      if(counts)
          outfilea <<(char*) this->BBNeighbors[i][1] << " " << betadist / (1.0 * (numneighbors))<< endl;
      else
          outfilea <<(char*) this->BBNeighbors[i][1] << " " << betadist / (1.0 * (numneighbors-1))<< endl;
  }
  outfilea.close();

  DestroyMetricSpace(MSQuery);
  fclose(AlphaFile);
  fclose(BetaFile);
  fclose(QueryFile);
  delete[] inputchar;
}
void lpGenerator::setshift(float s){
    setShift(s);
}

void lpGenerator::findNeighborsOneClassTwoDistanceAveLoad(char* AlphaIF, char* BetaIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand, bool counts) {
  const char* AlphaFilename = (const char* ) AlphaIF;
  const char* BetaFilename = (const char* ) BetaIF;
  const char* QueryFilename = (const char* ) QueryIF;
  string r = "r";
  FILE* AlphaFile = fopen(AlphaFilename, r.c_str());
  FILE* BetaFile = fopen(BetaFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());
  if (firstIteration == true){
      float epsilon = 1;
      float aarho = 0.9;
      float ssrho = 0.3;
      setEpsilonRho(epsilon,aarho, ssrho);
  }

  Distance adistfunc;
  Distance bdistfunc;
  cout << "building trees" << endl;
  if (firstIteration == true){
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
      nnfinder1.loadDispersionTree(AlphaFile, "trees/alpha" + to_string(rand) + "4.tree", "defaultaw",rand);
      nnfinder4.loadDispersionTree(BetaFile, "trees/beta" + to_string(rand) + "4.tree","defaultbw",rand);
      adistfunc = WeightedWordDist;
      bdistfunc = WeightedWordDist;
  }else{
      nnfinder1.loadDispersionTree(AlphaFile, "trees/alpha" + to_string(rand) + ".tree", "alpha",rand);
      nnfinder4.loadDispersionTree(BetaFile, "trees/beta" + to_string(rand) + ".tree","beta",rand);
      adistfunc = wordDist2Alpha;
      bdistfunc = wordDist2Beta;
  }
  /*this->buildalphatree( AlphaFile, 1, firstIteration, rand);
  this->buildbetatree(BetaFile, 1, firstIteration, rand);
  */
  cout << "I'm finding the neighbors and their two distances" << endl;

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  MetricSpace* MSQuery = CreateMetricSpace(inputchar, QueryFile);
  lpGenerator::TRSSize = MSQuery->numpoints;
  cout << "built the trees" << endl;

  nnfinder1.queryAlpha(numneighbors + 1, MSQuery, &this->AANeighbors, 1);
  nnfinder4.queryBeta(numneighbors + 1, MSQuery, &this->BBNeighbors, 1);
  cout << "all threads joined" << endl;
  string a = outputfile + "a";
  ofstream outfilea(a.c_str());
  for (int i = 0;i < TRSSize;i++) {
      float alphadist = 0;
      float betadist = 0;
      string inword( (char*) MSQuery->points[i]);
      if(counts){
          int anumsummed = 0;
          int aj =0;
          int acurcount = 0;
          int bnumsummed = 0;
          int bj =0;
          int bcurcount = 0;
        string q((char*)this->AANeighbors[i][aj]);
        if(this->countmap.find(q) != this->countmap.end()){
            acurcount = this->countmap[q];
        }
        else{
            acurcount = 1;
        }

        if(strcmp(inword.c_str(),(char*) this->AANeighbors[i][aj]) == 0)
            acurcount--;
        string p((char*)this->BBNeighbors[i][bj]);
        if(this->countmap.find(p) != this->countmap.end()){
            bcurcount = this->countmap[p];
        }
        else{
            bcurcount = 1;
        }

        if(strcmp(inword.c_str(),(char*) this->BBNeighbors[i][bj]) == 0)
            bcurcount--;
          while(anumsummed < numneighbors || bnumsummed < numneighbors){
              if (acurcount > 0 && anumsummed < numneighbors){
                  alphadist += adistfunc(MSQuery->points[i],this->AANeighbors[i][aj]);
                  acurcount = 0;
                  anumsummed++;
              }
              else{
                  aj++;
                  string u((char*)this->AANeighbors[i][aj]);
                  if(this->countmap.find(u) != this->countmap.end()){
                      acurcount = this->countmap[u];
                  }
                  else{
                      acurcount = 1;
                  }

                  if(strcmp(inword.c_str(),(char*) this->AANeighbors[i][aj]) == 0)
                      acurcount--;

              }
              if (bcurcount > 0 && bnumsummed < numneighbors){
                  betadist += bdistfunc(MSQuery->points[i],this->BBNeighbors[i][bj]);
                  bcurcount = 0;
                  bnumsummed++;
              }
              else{
                  bj++;
                  string t((char*)this->BBNeighbors[i][bj]);
                  if(this->countmap.find(t) != this->countmap.end()){
                      bcurcount = this->countmap[t];
                  }
                  else{
                      bcurcount = 1;
                  }
  
                  if(strcmp(inword.c_str(),(char*) this->BBNeighbors[i][bj]) == 0)
                      bcurcount--;
              }
          }

      }else{
          for (int j = 0; j < numneighbors - 1; j++){
              if (type == 1){
                  alphadist += adistfunc(MSQuery->points[i], this->AANeighbors[i][j+1]);
                  betadist += bdistfunc(MSQuery->points[i], this->BBNeighbors[i][j]);
                }
              else if (type == 2) {
                  alphadist += adistfunc(MSQuery->points[i], this->AANeighbors[i][j]);
                  betadist += bdistfunc(MSQuery->points[i], this->BBNeighbors[i][j+1]);
              }
              else if (type == 3) {
                  alphadist += adistfunc(MSQuery->points[i], this->AANeighbors[i][j]);
                  betadist += bdistfunc(MSQuery->points[i], this->BBNeighbors[i][j]);
              }
          }
      }
      outfilea << (char*) MSQuery->points[i] << " ";
      if(counts)
          outfilea <<(char*) this->AANeighbors[i][1] << " " << alphadist / (1.0 * (numneighbors)) << " " ;
      else
          outfilea <<(char*) this->AANeighbors[i][1] << " " << alphadist / (1.0 * (numneighbors-1)) << " " ;
      outfilea << (char*) MSQuery->points[i] << " ";
      if(counts)
          outfilea <<(char*) this->BBNeighbors[i][1] << " " << betadist / (1.0 * (numneighbors))<< endl;
      else
          outfilea <<(char*) this->BBNeighbors[i][1] << " " << betadist / (1.0 * (numneighbors-1))<< endl;
  }
  outfilea.close();

  DestroyMetricSpace(MSQuery);
  fclose(AlphaFile);
  fclose(BetaFile);
  fclose(QueryFile);
  delete[] inputchar;
}
void lpGenerator::findNeighborsOneClassTwoDistanceAveLoadCoil(char* AlphaIF, char* BetaIF,char* CoilIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand, bool counts) {
  const char* AlphaFilename = (const char* ) AlphaIF;
  const char* BetaFilename = (const char* ) BetaIF;
  const char* CoilFilename = (const char* ) CoilIF;
  const char* QueryFilename = (const char* ) QueryIF;
  string r = "r";
  FILE* AlphaFile = fopen(AlphaFilename, r.c_str());
  FILE* BetaFile = fopen(BetaFilename, r.c_str());
  FILE* CoilFile = fopen(CoilFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());
  if (firstIteration == true){
      float epsilon = 1;
      float aarho = 0.9;
      float ssrho = 0.3;
      setEpsilonRho(epsilon,aarho, ssrho);
  }

  Distance adistfunc;
  Distance bdistfunc;
  Distance cdistfunc;
  cout << "building trees" << endl;
  if (firstIteration == true){
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
      nnfinder1.loadDispersionTree(AlphaFile, "trees/alpha" + to_string(rand) + "4.tree", "defaultaw",rand);
      nnfinder4.loadDispersionTree(BetaFile, "trees/beta" + to_string(rand) + "4.tree","defaultbw",rand);
      nnfinder7.loadDispersionTree(CoilFile, "trees/coil" + to_string(rand) + "4.tree","defaultcw",rand);
      adistfunc = WeightedWordDist;
      bdistfunc = WeightedWordDist;
      cdistfunc = WeightedWordDist;
  }else{
      nnfinder1.loadDispersionTree(AlphaFile, "trees/alpha" + to_string(rand) + ".tree", "alpha",rand);
      nnfinder4.loadDispersionTree(BetaFile, "trees/beta" + to_string(rand) + ".tree","beta",rand);
      nnfinder7.loadDispersionTree(CoilFile, "trees/coil" + to_string(rand) + ".tree","coil",rand);
      adistfunc = wordDist2Alpha;
      bdistfunc = wordDist2Beta;
      cdistfunc = wordDist2Coil;
  }
  /*this->buildalphatree( AlphaFile, 1, firstIteration, rand);
  this->buildbetatree(BetaFile, 1, firstIteration, rand);
  */
  cout << "I'm finding the neighbors and their two distances" << endl;

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  MetricSpace* MSQuery = CreateMetricSpace(inputchar, QueryFile);
  lpGenerator::TRSSize = MSQuery->numpoints;
  cout << "built the trees" << endl;

  nnfinder1.queryAlpha(numneighbors + 1, MSQuery, &this->AANeighbors, 1);
  nnfinder4.queryBeta(numneighbors + 1, MSQuery, &this->BBNeighbors, 1);
  nnfinder7.queryCoil(numneighbors + 1, MSQuery, &this->CCNeighbors, 1);
  cout << "all threads joined" << endl;
  string a = outputfile + "a";
  ofstream outfilea(a.c_str());
  for (int i = 0;i < TRSSize;i++) {
      float alphadist = 0;
      float betadist = 0;
      float coildist = 0;
      string inword( (char*) MSQuery->points[i]);
      if(counts){
          int anumsummed = 0;
          int aj =0;
          int acurcount = 0;
          int bnumsummed = 0;
          int bj =0;
          int bcurcount = 0;
          int cnumsummed = 0;
          int cj =0;
          int ccurcount = 0;
        string q((char*)this->AANeighbors[i][aj]);
        if(this->countmap.find(q) != this->countmap.end()){
            acurcount = this->countmap[q];
        }
        else{
            acurcount = 1;
        }

        if(strcmp(inword.c_str(),(char*) this->AANeighbors[i][aj]) == 0)
            acurcount--;
        if (acurcount > 1)
            acurcount = 1;
        string p((char*)this->BBNeighbors[i][bj]);
        if(this->countmap.find(p) != this->countmap.end()){
            bcurcount = this->countmap[p];
        }
        else{
            bcurcount = 1;
        }
        if (bcurcount > 1)
            bcurcount = 1;

        if(strcmp(inword.c_str(),(char*) this->BBNeighbors[i][bj]) == 0)
            bcurcount--;
        string s((char*)this->CCNeighbors[i][cj]);
        if(this->countmap.find(s) != this->countmap.end()){
            ccurcount = this->countmap[s];
        }
        else{
            ccurcount = 1;
        }

        if(strcmp(inword.c_str(),(char*) this->CCNeighbors[i][cj]) == 0)
            ccurcount--;
        if (ccurcount > 1)
            ccurcount = 1;
        for (int d = 1; d <= numneighbors; d++){
      outfilea << (char*) MSQuery->points[i] << " ";
          while(anumsummed < d ){
              if (acurcount > 0 && anumsummed < d){
                  alphadist = adistfunc(MSQuery->points[i],this->AANeighbors[i][aj]);
          outfilea <<(char*) this->AANeighbors[i][aj] << " " << alphadist  << " " ;
                  acurcount = 0;
                  anumsummed++;
              }
              else{
                  aj++;
                  string u((char*)this->AANeighbors[i][aj]);
                  if(this->countmap.find(u) != this->countmap.end()){
                      acurcount = this->countmap[u];
                  }
                  else{
                      acurcount = 1;
                  }

                  if(strcmp(inword.c_str(),(char*) this->AANeighbors[i][aj]) == 0)
                      acurcount--;

              }
          }
          while(bnumsummed < d ){
              if (bcurcount > 0 && bnumsummed < d){
                  betadist = bdistfunc(MSQuery->points[i],this->BBNeighbors[i][bj]);
      outfilea << (char*) MSQuery->points[i] << " ";
          outfilea <<(char*) this->BBNeighbors[i][bj] << " " << betadist  << " " ;
                  bcurcount = 0;
                  bnumsummed++;
              }
              else{
                  bj++;
                  string t((char*)this->BBNeighbors[i][bj]);
                  if(this->countmap.find(t) != this->countmap.end()){
                      bcurcount = this->countmap[t];
                  }
                  else{
                      bcurcount = 1;
                  }
  
                  if(strcmp(inword.c_str(),(char*) this->BBNeighbors[i][bj]) == 0)
                      bcurcount--;
              }
          }
          while(cnumsummed < d ){
              if (ccurcount > 0 && cnumsummed < d){
                  coildist = cdistfunc(MSQuery->points[i],this->CCNeighbors[i][cj]);
      outfilea << (char*) MSQuery->points[i] << " ";
          outfilea <<(char*) this->CCNeighbors[i][cj] << " " << coildist  << endl ;
                  ccurcount = 0;
                  cnumsummed++;
              }
              else{
                  cj++;
                  string t((char*)this->CCNeighbors[i][cj]);
                  if(this->countmap.find(t) != this->countmap.end()){
                      ccurcount = this->countmap[t];
                  }
                  else{
                      ccurcount = 1;
                  }
  
                  if(strcmp(inword.c_str(),(char*) this->CCNeighbors[i][cj]) == 0)
                      ccurcount--;
              }
          }
        }

      }else{
          for (int j = 0; j < numneighbors - 1; j++){
              if (type == 1){
                  alphadist += adistfunc(MSQuery->points[i], this->AANeighbors[i][j+1]);
                  betadist += bdistfunc(MSQuery->points[i], this->BBNeighbors[i][j]);
                  coildist += cdistfunc(MSQuery->points[i], this->CCNeighbors[i][j]);
                }
              else if (type == 2) {
                  alphadist += adistfunc(MSQuery->points[i], this->AANeighbors[i][j]);
                  betadist += bdistfunc(MSQuery->points[i], this->BBNeighbors[i][j+1]);
                  coildist += cdistfunc(MSQuery->points[i], this->CCNeighbors[i][j]);
              }
              else if (type == 3) {
                  alphadist += adistfunc(MSQuery->points[i], this->AANeighbors[i][j]);
                  betadist += bdistfunc(MSQuery->points[i], this->BBNeighbors[i][j]);
                  coildist += cdistfunc(MSQuery->points[i], this->CCNeighbors[i][j+1]);
              }
          }
      }
      if(!counts){
      outfilea << (char*) MSQuery->points[i] << " ";
      if(counts)
          outfilea <<(char*) this->AANeighbors[i][1] << " " << alphadist / (1.0 * (numneighbors)) << " " ;
      else
          outfilea <<(char*) this->AANeighbors[i][1] << " " << alphadist / (1.0 * (numneighbors-1)) << " " ;
      outfilea << (char*) MSQuery->points[i] << " ";
      if(counts)
          outfilea <<(char*) this->BBNeighbors[i][1] << " " << betadist / (1.0 * (numneighbors))<< " ";
      else
          outfilea <<(char*) this->BBNeighbors[i][1] << " " << betadist / (1.0 * (numneighbors-1))<< " ";
      outfilea << (char*) MSQuery->points[i] << " ";
      if(counts)
          outfilea <<(char*) this->CCNeighbors[i][1] << " " << coildist / (1.0 * (numneighbors))<< endl;
      else
          outfilea <<(char*) this->CCNeighbors[i][1] << " " << coildist / (1.0 * (numneighbors-1))<< endl;
  }
  }
  outfilea.close();

  DestroyMetricSpace(MSQuery);
  fclose(AlphaFile);
  fclose(BetaFile);
  fclose(CoilFile);
  fclose(QueryFile);
  delete[] inputchar;
}

void lpGenerator::findNeighborsOneClassTwoDistanceWeighted(char* AlphaIF, char* BetaIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand, float rho) {
  const char* AlphaFilename = (const char* ) AlphaIF;
  const char* BetaFilename = (const char* ) BetaIF;
  const char* QueryFilename = (const char* ) QueryIF;
  string r = "r";
  FILE* AlphaFile = fopen(AlphaFilename, r.c_str());
  FILE* BetaFile = fopen(BetaFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());
  setEpsilonRho(1.0,rho, 0.3);
  string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
  readinweights(this->lengthOfWord,weightsfilename.c_str());

  cout << "building trees" << endl;
  this->buildalphatreeweighted( AlphaFile, 1, firstIteration, rand);
  this->buildbetatreeweighted(BetaFile, 1, firstIteration, rand);
  cout << "I'm finding the neighbors and their two distances" << endl;

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  MetricSpace* MSQuery = CreateMetricSpace(inputchar, QueryFile);
  lpGenerator::TRSSize = MSQuery->numpoints;
  cout << "built the trees" << endl;

  nnfinder1.queryAlpha(numneighbors, MSQuery, &this->AANeighbors, 1);
  nnfinder4.queryBeta(numneighbors, MSQuery, &this->BBNeighbors, 1);
  cout << "all threads joined" << endl;
  string a = outputfile + "a";
  ofstream outfilea(a.c_str());
  for (int i = 0;i < TRSSize;i++) {
      for (int j = 0; j < numneighbors - 1; j++){
          if (type == 1){
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->AANeighbors[i][j + 1] << " " << WeightedAlphaDistance(MSQuery->points[i], this->AANeighbors[i][j + 1]) << " " ;
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->BBNeighbors[i][j] << " " << WeightedBetaDistance(MSQuery->points[i], this->BBNeighbors[i][j]) << endl;
            }
          else if (type == 2) {
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->AANeighbors[i][j] << " " << WeightedAlphaDistance(MSQuery->points[i], this->AANeighbors[i][j]) << " " ;
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->BBNeighbors[i][j+1] << " " << WeightedBetaDistance(MSQuery->points[i], this->BBNeighbors[i][j+1]) << endl ;
          }
          else if (type == 3) {
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->AANeighbors[i][j] << " " << WeightedAlphaDistance(MSQuery->points[i], this->AANeighbors[i][j]) << " " ;
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->BBNeighbors[i][j] << " " << WeightedBetaDistance(MSQuery->points[i], this->BBNeighbors[i][j]) << endl;
          }
      }
  }
  outfilea.close();

  DestroyMetricSpace(MSQuery);
  fclose(AlphaFile);
  fclose(BetaFile);
  fclose(QueryFile);
  delete[] inputchar;

}

void lpGenerator::predict(char* AlphaIF, char* BetaIF, char* AlphaTS, char type, double alphatau, double betatau, string filename, int rand) {

  const char* AlphaInputFile = (const char* ) AlphaIF;
  const char* BetaInputFile = (const char* ) BetaIF;
  string r = "r";
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  bool firstIteration = false;

  thread alphatree(&lpGenerator::buildalphatree, this, alphaFile, 1, firstIteration,rand);
  thread betatree(&lpGenerator::buildbetatree, this, betaFile, 1, firstIteration,rand);
  alphatree.join();
  betatree.join();

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  FILE* AlphaTSFile = fopen(AlphaTS, r.c_str());
  MetricSpace* MSAlpha = CreateMetricSpace(inputchar, AlphaTSFile);
  lpGenerator::TRSSize = MSAlpha->numpoints;
  nnfinder1.queryAlpha(2, MSAlpha, &this->AANeighbors, 1);
  nnfinder4.queryBeta(2, MSAlpha, &this->BBNeighbors, 1);
  char predictedclass;
  ofstream outfile(filename.c_str());
  for (int i = 0; i < TRSSize; i++) {
	  double alphadist, betadist;
	  if (type == 'A') {
		   alphadist = wordDist2Alpha(MSAlpha->points[i], AANeighbors[i][1]);
		   betadist = wordDist2Beta(MSAlpha->points[i], BBNeighbors[i][0]);
           /*
		   if (alphadist < alphatau){
			   if (betadist < betatau){
				   if (alphadist < betadist){
					   predictedclass = 'A';
				   }
				   else {
					   predictedclass = 'B';
				   }
			   }
			   else {
				   predictedclass = 'A';
			   }
		   }
		   else {
			   if (betadist < betatau){
				predictedclass = 'B';
			   }
			   else {
				   predictedclass = 'C';
			   }
		   }
           */
		   /*
		  if(alphadist < betadist){
			//if (alphadist < alphatau)
			//	predictedclass = 'A';
			//else
			//	predictedclass = 'B';
			predictedclass = 'A';
			if (alphadist > coiltau)
				predictedclass = 'C';
			}
		else {
			predictedclass = 'B';
			if (betadist > coiltau)
				predictedclass = 'C';
			}
*/

		}
	  else if (type == 'B') {
		   alphadist = wordDist2Alpha(MSAlpha->points[i], AANeighbors[i][0]);
		   betadist = wordDist2Beta(MSAlpha->points[i], BBNeighbors[i][1]);
           /*BBNeighbors
		   if (alphadist < alphatau){
			   if (betadist < betatau){
				   if (alphadist < betadist){
					   predictedclass = 'A';
				   }
				   else {
					   predictedclass = 'B';
				   }
			   }
			   else {
				   predictedclass = 'A';
			   }
		   }
		   else {
			   if (betadist < betatau){
				predictedclass = 'B';
			   }
			   else {
				   predictedclass = 'C';
			   }
		   }
           */
		   /*
		  if(alphadist < betadist){
			//if (alphadist < alphatau)
			//	predictedclass = 'A';
			//else
			//	predictedclass = 'B';
			predictedclass = 'A';
			if (alphadist > coiltau)
				predictedclass = 'C';
			}
		else {
			predictedclass = 'B';
			if (betadist > coiltau)
				predictedclass = 'C';
			}
		*/

		}
	  else if (type == 'C') {
		   alphadist = wordDist2Alpha(MSAlpha->points[i], AANeighbors[i][0]);
		   betadist = wordDist2Beta(MSAlpha->points[i], BBNeighbors[i][0]);
           /*
		   if (alphadist < alphatau){
			   if (betadist < betatau){
				   if (alphadist < betadist){
					   predictedclass = 'A';
				   }
				   else {
					   predictedclass = 'B';
				   }
			   }
			   else {
				   predictedclass = 'A';
			   }
		   }
		   else {
			   if (betadist < betatau){
				predictedclass = 'B';
			   }
			   else {
				   predictedclass = 'C';
			   }
		   }
           */
		   /*
		  if(alphadist < betadist){
			//if (alphadist < alphatau)
			//	predictedclass = 'A';
			//else
			//	predictedclass = 'B';
			predictedclass = 'A';
			if (alphadist > coiltau)
				predictedclass = 'C';
			}
		else {
			predictedclass = 'B';
			if (betadist > coiltau)
				predictedclass = 'C';
			}

		*/
		}
           if (alphadist > alphatau && betadist > betatau)
               predictedclass = 'C';
           else{
               if (alphadist < betadist)
                   predictedclass = 'A';
               else
                   predictedclass = 'B';
           }
  	outfile << type << " " << predictedclass << " " << alphadist << " " << betadist << " " << endl;
	}

	outfile.close();
}
void lpGenerator::predictProteinLoadTree(char* AlphaIF, char* BetaIF, char* ProteinIF, double alphatau, double betatau, string filename, int numneighbors, int rand, string Treefile, string Treefile2, bool parallel, bool weights) {

  const char* AlphaInputFile = (const char* ) AlphaIF;
  const char* BetaInputFile = (const char* ) BetaIF;
  string r = "r";
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  cout << "###################################" << BetaIF << endl;
  if(weights == true){
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
      cout << "should have loaded in weights" << weightsfilename << this->lengthOfWord<< endl;
  }

  Distance alphadistfunc;
  Distance betadistfunc;
  if(weights == true){
      nnfinder1.loadDispersionTree(alphaFile, Treefile,"alphaweighted",rand);
      nnfinder4.loadDispersionTree(betaFile, Treefile2,"betaweighted",rand);
      alphadistfunc = WeightedAlphaDistance;
      betadistfunc = WeightedBetaDistance;
  }else{
      nnfinder1.loadDispersionTree(alphaFile, Treefile, "alpha",rand);
      nnfinder4.loadDispersionTree(betaFile, Treefile2,"beta",rand);
      alphadistfunc = wordDist2Alpha;
      betadistfunc = wordDist2Beta;
  }
  cout << "Trees loaded" << endl;
  /*thread alphatree(&lpGenerator::buildalphatree, this, alphaFile, 1, firstIteration,rand);
  thread betatree(&lpGenerator::buildbetatree, this, betaFile, 1, firstIteration,rand);
  alphatree.join();
  betatree.join();
  */

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  FILE* ProteinTSFile = fopen(ProteinIF, r.c_str());
  //at this point is when we need to split up the file and make each word into its own thing.  This is going to be the hardest part.
  const int num_splits = 8;
  int numprocessed = 0;
  MetricSpace* MSAlpha = CreateMetricSpace(inputchar, ProteinTSFile);
  if (parallel == true){
      thread* threads[num_splits];
      rewind(ProteinTSFile);
      MetricSpace** MSsAlpha = CreateMetricSpaces(inputchar, ProteinTSFile, num_splits);
      for (int i = 0; i < num_splits;i++){
          MetricSpace* MSAlpha2 = MSsAlpha[i];
          //threaddata td;
          //td.numneighbors = numneighbors;
          //td.MSAlpha = MSAlpha2;
          //td.lpgen = this;
          //td.treenum = 1;
          //td.writespace = this->AANeighbors + numprocessed;
          //pthread_create(&threads[i],NULL,(void*) lpGenerator::threadalphafunc, (void*)&td);

          thread t1(&lpGenerator::threadalphafunc, this, k + 1, MSAlpha, &this->AANeighbors + numprocessed, 1);
          threads[i] = &t1;
          numprocessed += MSAlpha2->numpoints;
      }
      for (int i =0;i < num_splits;i++){
          threads[i]->join();
      }


  }
  else{
      lpGenerator::TRSSize = MSAlpha->numpoints;
      //int numneighbors = 2;
      nnfinder1.queryAlpha(numneighbors, MSAlpha, &this->AANeighbors, 1);
      nnfinder4.queryBeta(numneighbors, MSAlpha, &this->BBNeighbors, 1);
  }
  char predictedclass;
  ofstream outfile(filename.c_str());
  for (int i = 0; i < TRSSize; i++) {
	  double alphadist = 0;
      double betadist = 0;
      for (int j = 0;j < numneighbors; j++){
		   alphadist += alphadistfunc(MSAlpha->points[i], AANeighbors[i][j]);
		   betadist += betadistfunc(MSAlpha->points[i], BBNeighbors[i][j]);
      }
      alphadist = alphadist / (float) numneighbors;
      betadist = betadist / (float) numneighbors;
		  if(alphadist <= betadist){
			predictedclass = 'A';
			if (alphadist > alphatau)
				predictedclass = 'C';
			}
		else {
			predictedclass = 'B';
			if (betadist > betatau)
				predictedclass = 'C';
			}

  	outfile << predictedclass << " " << alphadist << " " << betadist << " " << (char*)MSAlpha->points[i] << " " << (char*)AANeighbors[i][0] << " " << (char*)BBNeighbors[i][0] << endl;
	}

	outfile.close();
}

struct threaddata {
    int numneighbors;
    MetricSpace* MSAlpha;
    Pointer*** writespace;
    int treenum;
    lpGenerator* lpgen;
} ;

void* threadalphafunc2(void* tdp) {
    threaddata* td = (struct threaddata*) tdp;
    if (td->treenum == 1){
        cout << "before neighbor finding" << endl;
        td->lpgen->nnfinder1.queryAlpha(td->numneighbors, td->MSAlpha, td->writespace, td->treenum);
    }
    else if (td->treenum == 2)
        td->lpgen->nnfinder2.queryAlpha(td->numneighbors, td->MSAlpha, td->writespace, td->treenum);
    else if (td->treenum == 3)
        td->lpgen->nnfinder3.queryAlpha(td->numneighbors, td->MSAlpha, td->writespace, td->treenum);
    cout << "donethreadalphafunc" << endl;
    return NULL;
}

void* threadbetafunc2(void* tdp) {
    threaddata* td = (struct threaddata*) tdp;
    if (td->treenum == 1)
        td->lpgen->nnfinder4.queryBeta(td->numneighbors, td->MSAlpha, td->writespace, td->treenum);
    else if (td->treenum == 2)
        td->lpgen->nnfinder5.queryBeta(td->numneighbors, td->MSAlpha, td->writespace, td->treenum);
    else if (td->treenum == 3)
        td->lpgen->nnfinder6.queryBeta(td->numneighbors, td->MSAlpha, td->writespace, td->treenum);
    cout << "donethreadalphafunc" << endl;
    return NULL;
}

void lpGenerator::predictProteinLoadTreeDefault(char* AlphaIF, char* BetaIF, char* ProteinIF, double alphatau, double betatau, string filename, int numneighbors, int rand, string Treefile1,string Treefile2, bool parallel, bool weights,bool training) {

  const char* AlphaInputFile = (const char* ) AlphaIF;
  const char* BetaInputFile = (const char* ) BetaIF;
  string r = "r";
  setEpsilonRho(1,0.9,0.3);
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  cout << "###################################" << BetaIF << endl;

  Distance alphadistfunc;
  Distance betadistfunc;
  if(weights == true){
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
      cout << "should have loaded in weights" << weightsfilename << this->lengthOfWord<< endl;
  }
  if(weights){
      nnfinder1.loadDispersionTree(alphaFile, Treefile1, "defaultaw",rand);
      nnfinder4.loadDispersionTree(betaFile, Treefile2,"defaultbw",rand);
      alphadistfunc = WeightedWordDist;
      betadistfunc = WeightedWordDist;
  }else{
  nnfinder1.loadDispersionTree(alphaFile, Treefile1, "defaulta",rand);
  nnfinder4.loadDispersionTree(betaFile, Treefile2,"defaultb",rand);
      alphadistfunc = WordDist;
      betadistfunc = WordDist;
  }
  cout << "Trees loaded" << endl;
  /*thread alphatree(&lpGenerator::buildalphatree, this, alphaFile, 1, firstIteration,rand);
  thread betatree(&lpGenerator::buildbetatree, this, betaFile, 1, firstIteration,rand);
  alphatree.join();
  betatree.join();
  */

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  FILE* ProteinTSFile = fopen(ProteinIF, r.c_str());
  //at this point is when we need to split up the file and make each word into its own thing.  This is going to be the hardest part.
  const int num_splits = 2;
  int numprocessed = 0;
  MetricSpace* MSAlpha = CreateMetricSpace(inputchar, ProteinTSFile);
  if (parallel == true){
      pthread_t threads[num_splits];
      pthread_t bthreads[num_splits];
      rewind(ProteinTSFile);
      MetricSpace** MSsAlpha = CreateMetricSpaces(inputchar, ProteinTSFile, num_splits);
      threaddata* tds = (threaddata*)  malloc(num_splits * sizeof(threaddata));
      for (int i = 0; i < num_splits;i++){
          MetricSpace* MSAlpha2 = MSsAlpha[i];
          threaddata td;
          td.numneighbors = numneighbors;
          td.MSAlpha = MSsAlpha[i];
          td.lpgen = this;
          td.treenum = 1;
          td.writespace = &this->AANeighbors + numprocessed;
          tds[i] = td;
          //threaddata td2;
          //td2.numneighbors = numneighbors;
          //td2.MSAlpha = MSAlpha2;
          //td2.lpgen = this;
          //td2.treenum = 1;
          //td2.writespace = &this->BBNeighbors + numprocessed;
          pthread_create(&threads[i],NULL,threadalphafunc2, (void*)&tds[i]);
          //pthread_create(&bthreads[i],NULL,threadbetafunc2, (void*)&td2);

          //thread t1(&lpGenerator::threadalphafunc, this, numneighbors, MSsAlpha[i], &this->AANeighbors + numprocessed, 1);
          //thread t2(&lpGenerator::threadbetafunc, this, numneighbors, MSsAlpha[i], &this->BBNeighbors + numprocessed, 1);
          cout << "spawned threads" << endl;
          //threads[i] = &t1;
          //bthreads[i] = &t2;
          cout << "added the theads" << endl;
          numprocessed += MSAlpha2->numpoints;
          cout << "processed the stuff " << td.MSAlpha->numpoints << endl;
      }
      cout << "got out of the spawning" << endl;
      for (int i =0;i < num_splits;i++){
          pthread_join(threads[i],NULL);
          cout << "threads joined" << endl;
          //pthread_join(bthreads[i],NULL);
          //threads[i]->join();
      }
      numprocessed = 0;
      
      threaddata* tds2 = (threaddata*)  malloc(num_splits * sizeof(threaddata));
      for (int i = 0; i < num_splits;i++){
          MetricSpace* MSAlpha2 = MSsAlpha[i];
          threaddata td2;
          td2.numneighbors = numneighbors;
          td2.MSAlpha = MSsAlpha[i];
          td2.lpgen = this;
          td2.treenum = 1;
          td2.writespace = &this->BBNeighbors + numprocessed;
          tds2[i] = td2;
          pthread_create(&bthreads[i],NULL,threadbetafunc2, (void*)&tds2[i]);

          numprocessed += MSAlpha2->numpoints;
      }
      for (int i =0;i < num_splits;i++){
          //pthread_join(threads[i],NULL);
          pthread_join(bthreads[i],NULL);
          cout << "threads b joined" << endl;

      }
      


  }
  else{
      //int numneighbors = 2;
      nnfinder1.queryAlpha(numneighbors, MSAlpha, &this->AANeighbors, 1);
      nnfinder4.queryBeta(numneighbors, MSAlpha, &this->BBNeighbors, 1);
  }
  lpGenerator::TRSSize = MSAlpha->numpoints;
  cout << "all stuff joined" << endl;
  char predictedclass;
  ofstream outfile(filename.c_str());
  for (int i = 0; i < TRSSize; i++) {
	  double alphadist = 0;
      double betadist = 0;
      for (int j = 0;j < numneighbors; j++){
          if (training){
            outfile << alphadistfunc(MSAlpha->points[i], AANeighbors[i][j]) << " " << betadistfunc(MSAlpha->points[i], BBNeighbors[i][j]) << " " << (char*)MSAlpha->points[i] << " " << (char*)AANeighbors[i][j] << " " << (char*)BBNeighbors[i][j] << endl;
          }
		   alphadist += alphadistfunc(MSAlpha->points[i], AANeighbors[i][j]);
		   betadist += betadistfunc(MSAlpha->points[i], BBNeighbors[i][j]);
      }
      alphadist = alphadist / (float) numneighbors;
      betadist = betadist / (float) numneighbors;
		  if(alphadist <= betadist){
			predictedclass = 'A';
			if (alphadist > alphatau)
				predictedclass = 'C';
			}
		else {
			predictedclass = 'B';
			if (betadist > betatau)
				predictedclass = 'C';
			}

        if(!training)
            outfile << predictedclass << " " << alphadist << " " << betadist << " " << (char*)MSAlpha->points[i] << " " << (char*)AANeighbors[i][0] << " " << (char*)BBNeighbors[i][0] << endl;
	}

	outfile.close();
}
void lpGenerator::predictProteinLoadTreeDefaultCoil(char* AlphaIF, char* BetaIF,char* CoilIF, char* ProteinIF, double alphatau, double betatau, string filename, int numneighbors, int rand, string Treefile1,string Treefile2,string Treefile3, bool parallel, bool weights,bool training) {

  const char* AlphaInputFile = (const char* ) AlphaIF;
  const char* BetaInputFile = (const char* ) BetaIF;
  const char* CoilInputFile = (const char* ) CoilIF;
  string r = "r";
  setEpsilonRho(1,0.9,0.3);
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  FILE* coilFile = fopen(CoilInputFile, r.c_str());
  cout << "###################################" << BetaIF << endl;

  Distance alphadistfunc;
  Distance betadistfunc;
  Distance coildistfunc;
  if(weights == true){
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
      cout << "should have loaded in weights" << weightsfilename << this->lengthOfWord<< endl;
  }
  if(weights){
      nnfinder1.loadDispersionTree(alphaFile, Treefile1, "defaultaw",rand);
      nnfinder4.loadDispersionTree(betaFile, Treefile2,"defaultbw",rand);
      nnfinder7.loadDispersionTree(coilFile, Treefile3,"defaultcw",rand);
      alphadistfunc = WeightedWordDist;
      betadistfunc = WeightedWordDist;
      coildistfunc = WeightedWordDist;
  }else{
  nnfinder1.loadDispersionTree(alphaFile, Treefile1, "defaulta",rand);
  nnfinder4.loadDispersionTree(betaFile, Treefile2,"defaultb",rand);
      alphadistfunc = WordDist;
      betadistfunc = WordDist;
  }
  cout << "Trees loaded" << endl;
  /*thread alphatree(&lpGenerator::buildalphatree, this, alphaFile, 1, firstIteration,rand);
  thread betatree(&lpGenerator::buildbetatree, this, betaFile, 1, firstIteration,rand);
  alphatree.join();
  betatree.join();
  */

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  FILE* ProteinTSFile = fopen(ProteinIF, r.c_str());
  //at this point is when we need to split up the file and make each word into its own thing.  This is going to be the hardest part.
  MetricSpace* MSAlpha = CreateMetricSpace(inputchar, ProteinTSFile);
  //int numneighbors = 2;
  nnfinder1.queryAlpha(numneighbors, MSAlpha, &this->AANeighbors, 1);
  nnfinder4.queryBeta(numneighbors, MSAlpha, &this->BBNeighbors, 1);
  nnfinder7.queryCoil(numneighbors, MSAlpha, &this->CCNeighbors, 1);
  lpGenerator::TRSSize = MSAlpha->numpoints;
  cout << "all stuff joined" << endl;
  char predictedclass;
  ofstream outfile(filename.c_str());
  for (int i = 0; i < TRSSize; i++) {
	  double alphadist = 0;
      double betadist = 0;
      double coildist = 0;
      for (int j = 0;j < numneighbors; j++){
          if (training){
            outfile << alphadistfunc(MSAlpha->points[i], AANeighbors[i][j]) << " " << betadistfunc(MSAlpha->points[i], BBNeighbors[i][j]) << " " << coildistfunc(MSAlpha->points[i],CCNeighbors[i][j]) << " " << (char*)MSAlpha->points[i] << " " << (char*)AANeighbors[i][j] << " " << (char*)BBNeighbors[i][j] << " " << (char*)CCNeighbors[i][j] << endl;
          }
		   alphadist += alphadistfunc(MSAlpha->points[i], AANeighbors[i][j]);
		   betadist += betadistfunc(MSAlpha->points[i], BBNeighbors[i][j]);
		   coildist += coildistfunc(MSAlpha->points[i], CCNeighbors[i][j]);
      }
      alphadist = alphadist / (float) numneighbors;
      betadist = betadist / (float) numneighbors;
      coildist = coildist / (float) numneighbors;
		  if(alphadist <= betadist && alphadist <= coildist){
			predictedclass = 'A';
			if (alphadist > alphatau)
				predictedclass = 'C';
			}
		else if (betadist <= coildist){
			predictedclass = 'B';
			if (betadist > betatau)
				predictedclass = 'C';
			}
        else
            predictedclass = 'C';

        if(!training)
            outfile << predictedclass << " " << alphadist << " " << betadist << " " << coildist << " "<< (char*)MSAlpha->points[i] << " " << (char*)AANeighbors[i][0] << " " << (char*)BBNeighbors[i][0] << " " << (char*)CCNeighbors[i][0] << endl;
	}

	outfile.close();
}
void lpGenerator::predictProteinLoadTreeCoil(char* AlphaIF, char* BetaIF,char* CoilIF, char* ProteinIF, double alphatau, double betatau, string filename, int numneighbors, int rand, string Treefile1,string Treefile2,string Treefile3, bool parallel, bool weights,bool training) {

  const char* AlphaInputFile = (const char* ) AlphaIF;
  const char* BetaInputFile = (const char* ) BetaIF;
  const char* CoilInputFile = (const char* ) CoilIF;
  string r = "r";
  setEpsilonRho(1,0.9,0.3);
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  FILE* coilFile = fopen(CoilInputFile, r.c_str());
  cout << "###################################" << BetaIF << endl;

  Distance alphadistfunc;
  Distance betadistfunc;
  Distance coildistfunc;
  if(weights == true){
      string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
      readinweights(this->lengthOfWord,weightsfilename.c_str());
      cout << "should have loaded in weights" << weightsfilename << this->lengthOfWord<< endl;
  }
  if(weights){
      nnfinder1.loadDispersionTree(alphaFile, Treefile1, "alphaweighted",rand);
      nnfinder4.loadDispersionTree(betaFile, Treefile2,"betaweighted",rand);
      nnfinder7.loadDispersionTree(coilFile, Treefile3,"coilweighted",rand);
      alphadistfunc = WeightedAlphaDistance;
      betadistfunc = WeightedBetaDistance;
      coildistfunc = WeightedCoilDistance;
  }else{
  nnfinder1.loadDispersionTree(alphaFile, Treefile1, "alpha",rand);
  nnfinder4.loadDispersionTree(betaFile, Treefile2,"beta",rand);
  nnfinder7.loadDispersionTree(coilFile, Treefile3,"coil",rand);
      alphadistfunc = wordDist2Alpha;
      betadistfunc = wordDist2Beta;
      coildistfunc = wordDist2Coil;
  }
  cout << "Trees loaded" << endl;
  /*thread alphatree(&lpGenerator::buildalphatree, this, alphaFile, 1, firstIteration,rand);
  thread betatree(&lpGenerator::buildbetatree, this, betaFile, 1, firstIteration,rand);
  alphatree.join();
  betatree.join();
  */

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  FILE* ProteinTSFile = fopen(ProteinIF, r.c_str());
  //at this point is when we need to split up the file and make each word into its own thing.  This is going to be the hardest part.
  MetricSpace* MSAlpha = CreateMetricSpace(inputchar, ProteinTSFile);
  //int numneighbors = 2;
  nnfinder1.queryAlpha(numneighbors, MSAlpha, &this->AANeighbors, 1);
  nnfinder4.queryBeta(numneighbors, MSAlpha, &this->BBNeighbors, 1);
  nnfinder7.queryCoil(numneighbors, MSAlpha, &this->CCNeighbors, 1);
  lpGenerator::TRSSize = MSAlpha->numpoints;
  cout << "all stuff joined" << endl;
  char predictedclass;
  ofstream outfile(filename.c_str());
  for (int i = 0; i < TRSSize; i++) {
	  double alphadist = 0;
      double betadist = 0;
      double coildist = 0;
      for (int j = 0;j < numneighbors; j++){
          if (training){
            outfile << alphadistfunc(MSAlpha->points[i], AANeighbors[i][j]) << " " << betadistfunc(MSAlpha->points[i], BBNeighbors[i][j]) << " " << coildistfunc(MSAlpha->points[i],CCNeighbors[i][j]) << " " << (char*)MSAlpha->points[i] << " " << (char*)AANeighbors[i][j] << " " << (char*)BBNeighbors[i][j] << " " << (char*)CCNeighbors[i][j] << endl;
          }
		   alphadist += alphadistfunc(MSAlpha->points[i], AANeighbors[i][j]);
		   betadist += betadistfunc(MSAlpha->points[i], BBNeighbors[i][j]);
		   coildist += coildistfunc(MSAlpha->points[i], CCNeighbors[i][j]);
      }
      alphadist = alphadist / (float) numneighbors;
      betadist = betadist / (float) numneighbors;
      coildist = coildist / (float) numneighbors;
		  if(alphadist <= betadist && alphadist <= coildist){
			predictedclass = 'A';
			if (alphadist > alphatau)
				predictedclass = 'C';
			}
		else if (betadist <= coildist){
			predictedclass = 'B';
			if (betadist > betatau)
				predictedclass = 'C';
			}
        else
            predictedclass = 'C';

        if(!training)
            outfile << predictedclass << " " << alphadist << " " << betadist << " " << coildist << " "<< (char*)MSAlpha->points[i] << " " << (char*)AANeighbors[i][0] << " " << (char*)BBNeighbors[i][0] << " " << (char*)CCNeighbors[i][0] << endl;
	}

	outfile.close();
}
void lpGenerator::predictProtein(char* AlphaIF, char* BetaIF, char* ProteinIF, double alphatau, double betatau, string filename, int numneighbors, int rand, bool firstIteration) {

  const char* AlphaInputFile = (const char* ) AlphaIF;
  const char* BetaInputFile = (const char* ) BetaIF;
  string r = "r";
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  setEpsilonRho(1,0.9,0.3);
  cout << "###################################" << BetaIF << endl;

  thread alphatree(&lpGenerator::buildalphatree, this, alphaFile, 1, firstIteration,rand);
  thread betatree(&lpGenerator::buildbetatree, this, betaFile, 1, firstIteration,rand);
  alphatree.join();
  betatree.join();

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  FILE* ProteinTSFile = fopen(ProteinIF, r.c_str());
  //at this point is when we need to split up the file and make each word into its own thing.  This is going to be the hardest part.
  MetricSpace* MSAlpha = CreateMetricSpace(inputchar, ProteinTSFile);
  lpGenerator::TRSSize = MSAlpha->numpoints;
  //int numneighbors = 2;
  nnfinder1.queryAlpha(numneighbors, MSAlpha, &this->AANeighbors, 1);
  nnfinder4.queryBeta(numneighbors, MSAlpha, &this->BBNeighbors, 1);
  char predictedclass;
  ofstream outfile(filename.c_str());
  for (int i = 0; i < TRSSize; i++) {
	  double alphadist = 0;
      double betadist = 0;
      for (int j = 0;j < numneighbors; j++){
           if (firstIteration == false){
               alphadist += wordDist2Alpha(MSAlpha->points[i], AANeighbors[i][j]);
               betadist += wordDist2Beta(MSAlpha->points[i], BBNeighbors[i][j]);
           }
           else{
               alphadist += WordDist(MSAlpha->points[i], AANeighbors[i][j]);
               betadist += WordDist(MSAlpha->points[i], BBNeighbors[i][j]);
           }
               
      }
      alphadist = alphadist / (float) numneighbors;
      betadist = betadist / (float) numneighbors;
		  if(alphadist <= betadist){
			predictedclass = 'A';
			if (alphadist > alphatau)
				predictedclass = 'C';
			}
		else {
			predictedclass = 'B';
			if (betadist > betatau)
				predictedclass = 'C';
			}

  	outfile << predictedclass << " " << alphadist << " " << betadist << " " << endl;
	}

	outfile.close();
}
void lpGenerator::predictProteinWithKNN(char* AlphaIF, char* BetaIF, char* ProteinIF, double alphatau, double betatau, string filename, int numneighbors, int rand) {

  const char* AlphaInputFile = (const char* ) AlphaIF;
  const char* BetaInputFile = (const char* ) BetaIF;
  string r = "r";
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  bool firstIteration = false;
  cout << "###################################" << BetaIF << endl;

  thread alphatree(&lpGenerator::buildalphatree, this, alphaFile, 1, firstIteration,rand);
  thread betatree(&lpGenerator::buildbetatree, this, betaFile, 1, firstIteration,rand);
  alphatree.join();
  betatree.join();

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  FILE* ProteinTSFile = fopen(ProteinIF, r.c_str());
  //at this point is when we need to split up the file and make each word into its own thing.  This is going to be the hardest part.
  MetricSpace* MSAlpha = CreateMetricSpace(inputchar, ProteinTSFile);
  lpGenerator::TRSSize = MSAlpha->numpoints;
  //int numneighbors = 2;
  nnfinder1.queryAlpha(numneighbors, MSAlpha, &this->AANeighbors, 1);
  nnfinder4.queryBeta(numneighbors, MSAlpha, &this->BBNeighbors, 1);
  ofstream outfile(filename.c_str());
  float alphadist;
  float betadist;
  for (int i = 0; i < TRSSize; i++) {
      for (int j = 0;j < numneighbors; j++){
		   alphadist = wordDist2Alpha(MSAlpha->points[i], AANeighbors[i][j]);
		   betadist = wordDist2Beta(MSAlpha->points[i], BBNeighbors[i][j]);
           outfile << (char*)MSAlpha->points[i] << " " << (char*)AANeighbors[i][j] << " " << alphadist << " " << (char*)BBNeighbors[i][j] << " " << betadist << endl;
      }

	}

	outfile.close();
}
void lpGenerator::predictProteinWithKNNCoil(char* AlphaIF, char* BetaIF,char* CoilIF, char* ProteinIF, string filename, int numneighbors, int rand, bool firstIteration, int dirnum,bool ss,bool weights,bool dirr) {
  printf("starting predictproteinwightknncoil with ss %d, weights %d, and dirr %d",ss,weights,dirr);

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
          cout << "using SecondaryDist" << endl;
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
          printf("Destroyed the beta tree\n");
          coilnn = nnfinder7.queryCoilRange(treesize, MSAlpha, &this->CCNeighbors, 1);
          printf("before Destroyed the coil tree\n");
          nnfinder7.destroyDTree('C');
          printf("Destroyed the coil tree\n");
         /* 
          nnfinder1.loadDispersionTree(alphaFile, "a.tree", "ass",rand);
          nnfinder4.loadDispersionTree(betaFile, "b.tree","bss",rand);
          nnfinder7.loadDispersionTree(coilFile, "c.tree","css",rand);
          */
      } else {
          alphadistfunc = WeightedWordDist2;
          betadistfunc = WeightedWordDist2;
          coildistfunc = WeightedWordDist2;
          nnfinder1.loadDispersionTree(alphaFile, "trees/alpha" + to_string(rand) + to_string(dirnum) + ".tree", "defaultaw",rand);
          nnfinder4.loadDispersionTree(betaFile, "trees/beta" + to_string(rand) + to_string(dirnum) + ".tree","defaultbw",rand);
          nnfinder7.loadDispersionTree(coilFile, "trees/coil" + to_string(rand) + to_string(dirnum) + ".tree","defaultcw",rand);
          nnfinder1.queryAlpha(numneighbors, MSAlpha, &this->AANeighbors, 1);
          nnfinder1.destroyDTree('A');
          nnfinder4.queryBeta(numneighbors, MSAlpha, &this->BBNeighbors, 1);
          nnfinder4.destroyDTree('B');
          printf("Destroyed the beta tree\n");
          nnfinder7.queryCoil(numneighbors, MSAlpha, &this->CCNeighbors, 1);
          printf("before Destroyed the coil tree\n");
          nnfinder7.destroyDTree('C');
          printf("Destroyed the coil tree\n");
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
          nnfinder1.loadDispersionTree(alphaFile, "trees/alpha" + to_string(rand) + ".tree", "alpha",rand);
          nnfinder4.loadDispersionTree(betaFile, "trees/beta" + to_string(rand) + ".tree","beta",rand);
          nnfinder7.loadDispersionTree(coilFile, "trees/coil" + to_string(rand) + ".tree","coil",rand);
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
      printf("alphann: %d betann: %d coilnn: %d",alphann,betann,coilnn);
      printf("%s\n",AANeighbors[0][0]);
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
      printf("finished the alphas\n");
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



void lpGenerator::findNeighbors(char* AlphaIF, char* BetaIF, char* AlphaTS, char* BetaTS, char* CoilTS, bool firstIteration, int rand) {
  const char* AlphaInputFile = (const char* ) AlphaIF;
  const char* BetaInputFile = (const char* ) BetaIF;
  string r = "r";
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  FILE* alphaFile2 = fopen(AlphaInputFile, r.c_str());
  FILE* alphaFile3 = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile2 = fopen(BetaInputFile, r.c_str());
  FILE* betaFile3 = fopen(BetaInputFile, r.c_str());

  thread alphatree(&lpGenerator::buildalphatree, this, alphaFile, 1, firstIteration,rand);
  thread betatree(&lpGenerator::buildbetatree, this, betaFile, 1, firstIteration,rand);
  thread alphatree2(&lpGenerator::buildalphatree, this, alphaFile2, 2, firstIteration,rand);
  thread betatree2(&lpGenerator::buildbetatree, this, betaFile2, 2, firstIteration,rand);
  thread alphatree3(&lpGenerator::buildalphatree, this, alphaFile3, 3, firstIteration,rand);
  thread betatree3(&lpGenerator::buildbetatree, this, betaFile3, 3, firstIteration,rand);
  alphatree3.join();
  betatree3.join();
  alphatree2.join();
  betatree2.join();
  alphatree.join();
  betatree.join();

  string input = "PROTEIN_WORDS";
  char* inputchar = new char[input.size() + 1];
  copy(input.begin(), input.end(), inputchar);
  inputchar[input.size()] = '\0';

  FILE* AlphaTSFile = fopen(AlphaTS, r.c_str());
  FILE* BetaTSFile = fopen(BetaTS, r.c_str());
  FILE* CoilTSFile = fopen(CoilTS, r.c_str());
  MetricSpace *MSBeta = CreateMetricSpace(inputchar, BetaTSFile);
  MetricSpace* MSAlpha = CreateMetricSpace(inputchar, AlphaTSFile);
  MetricSpace *MSCoil = CreateMetricSpace(inputchar, CoilTSFile);
  rewind(AlphaTSFile);
  rewind(BetaTSFile);
  rewind(CoilTSFile);
  MetricSpace *MSBeta2 = CreateMetricSpace(inputchar, BetaTSFile);
  MetricSpace* MSAlpha2 = CreateMetricSpace(inputchar, AlphaTSFile);
  MetricSpace *MSCoil2 = CreateMetricSpace(inputchar, CoilTSFile);
  lpGenerator::TRSSize = MSAlpha->numpoints;
  cout << "built the trees" << endl;

  //thread t1(&lpGenerator::threadalphafunc, this, k + 1, MSAlpha, &this->AANeighbors, 1);
  //thread t2(&lpGenerator::threadalphafunc, this,l, MSBeta2, &this->ABNeighbors, 2);
  //thread t3(&lpGenerator::threadbetafunc, this,k + 1, MSBeta, &this->BBNeighbors, 1);
  //thread t4(&lpGenerator::threadbetafunc, this, l, MSAlpha2, &this->BANeighbors, 2);
  //thread t5(&lpGenerator::threadalphafunc, this,l, MSCoil, &this->ACNeighbors, 3);
  //thread t6(&lpGenerator::threadbetafunc, this,l, MSCoil2, &this->BCNeighbors, 3);

  //t1.join();
  //t2.join();
  //t3.join();
  //t4.join();
  //t5.join();
  //t6.join();
  nnfinder1.queryAlpha(k+1, MSAlpha, &this->AANeighbors, 1);
  nnfinder2.queryAlpha(l, MSBeta, &this->BANeighbors, 2);
  nnfinder3.queryAlpha(l, MSCoil, &this->ACNeighbors, 3);
  nnfinder4.queryBeta(k+1, MSBeta2, &this->BBNeighbors, 1);
  nnfinder5.queryBeta(l, MSAlpha2, &this->ABNeighbors, 2);
  nnfinder6.queryBeta(l, MSCoil2, &this->BCNeighbors, 3);
  cout << "all threads joined" << endl;
  cout << (char*) AANeighbors[0][0] << " " << (char* ) AANeighbors[0][1] << " " << (char* ) AANeighbors[0][2] << endl;
  cout << "the distances are: " << WordDist(AANeighbors[0][0], AANeighbors[0][1]) << " and " << WordDist(AANeighbors[0][0], AANeighbors[0][2]) << endl;

  DestroyMetricSpace(MSAlpha);
  DestroyMetricSpace(MSBeta);
  DestroyMetricSpace(MSCoil);
  DestroyMetricSpace(MSAlpha2);
  DestroyMetricSpace(MSBeta2);
  DestroyMetricSpace(MSCoil2);
  fclose(AlphaTSFile);
  fclose(BetaTSFile);
  fclose(CoilTSFile);
  fclose(alphaFile);
  fclose(betaFile);
  delete[] inputchar;
}

void lpGenerator::createMapping() {
  Mapping m;
  this->mapping = &m;
  int VarNum = numClasses * lengthOfWord * (ALPHABET_SIZE * (ALPHABET_SIZE + 1) / 2);
  m.indexToClass = new int[VarNum + 1];
  m.indexToPosition = new int[VarNum + 1];
  m.indexToFirstLetter = new int[VarNum + 1];
  m.indexToSecondLetter = new int[VarNum + 1];
  m.variableToIndex = new int***[numClasses + 1];


  for( int i = 1; i <= numClasses; i++) {
    m.variableToIndex[i] = new int**[lengthOfWord + 1];
  }

  for(int i = 1; i <= numClasses; i++) {
    for(int j = 1;j <= lengthOfWord; j++) {
      m.variableToIndex[i][j] = new int*[ALPHABET_SIZE + 1];
    }
  }

  for(int i = 1; i <= numClasses; i++) {
    for(int j = 1;j <= lengthOfWord; j++) {
      for(int k = 1; k <= ALPHABET_SIZE; k++) {
        m.variableToIndex[i][j][k] = new int[ALPHABET_SIZE + 1];
      }
    }
  }

  int index = 0;
  for (int i = 1; i <= numClasses; i++) {
     for (int j = 1; j <= lengthOfWord; j++) {
        for (int k = 1; k <= ALPHABET_SIZE; k++) {
           for(int l = k; l <= ALPHABET_SIZE; l++) {
              m.variableToIndex[i][j][k][l] = ++index;
              m.variableToIndex[i][j][l][k] = index;    // residues2 at the same position are symmetric
              m.indexToClass[index] = i;
              m.indexToPosition[index] = j;
              m.indexToFirstLetter[index] = k;
              m.indexToSecondLetter[index] = l;
           }
        }
     }
  }

}

void lpGenerator::writeObjectiveFunction(string alphaTSname, string betaTSname, string coilTSname, ofstream* output, int t_ssize) {
  double alp_ab = 1.0/3;
  cout << t_ssize << endl;
  int a_querynum = t_ssize;
  int b_querynum = t_ssize;
  int c_querynum = t_ssize;
  double a_coef = (1.0/k) * (1.0/a_querynum) * alp_ab;
  double a_coef_imposters = (1.0/a_querynum) * (alp_ab) ;
  double b_coef = (1.0/k) * (1.0/b_querynum) * alp_ab ;
  double b_coef_imposters = (1.0/b_querynum) * (alp_ab) ;
  double c_coef = (1.0/c_querynum) * (1.0 / 2);
  ifstream awords(alphaTSname.c_str());
  ifstream bwords(betaTSname.c_str());
  ifstream cwords(coilTSname.c_str());
  string queryWord;
  string queryWordB;
  for(int i = 0; i < 3; i++) {
	  awords >> queryWord;
	  bwords >> queryWord;
	  cwords >> queryWord;
  }
  (*(output)) << "Minimize\n" << "obj: ";
  int i;
  int j;
  for (i = 0; i < a_querynum; i++) {
     awords >> queryWord;

     for (j = 0; j < k; j++) {
        if (i == 0 && j == 0)
          (*(output)) << " " <<   scientific << a_coef << " d_" << queryWord << "_" << j;
        else
          (*(output)) << " + "   << scientific << a_coef << " d_" << queryWord << "_" << j;
     }
     (*(output)) << " + " <<   scientific << a_coef_imposters << " d_" << queryWord;
     for (j = 0;j < k; j++) {
	     (*(output)) << " + " << scientific << a_coef << " d_" << queryWord << "_F" << j;
     }

     if (i%9 == 0) {
        (*(output)) << endl;
     }
  }
  for (i = 0; i < b_querynum; i++) {
     bwords >> queryWord;

     for (j = 0; j < k; j++) {
       (*(output)) << " + "   << scientific << b_coef << " d_" << queryWord << "_" << j;
     }
     (*(output)) << " + "   << scientific << b_coef_imposters << " d_" << queryWord;
     for (j = 0; j < k; j++) {
	     (*(output)) << " + " << scientific << b_coef << " d_" << queryWord << "_F" << j;
     }

     if (i%9 == 0) {
        (*(output)) << endl;
     }
  }

  for (i = 0; i < c_querynum; i++) {
     cwords >> queryWord;

     (*(output)) << " + "   << scientific << c_coef << " d_" << queryWord << "_A";
     (*(output)) << " + "   << scientific << c_coef << " d_" << queryWord << "_B";

     if (i%9 == 0) {
        (*(output)) << endl;
     }
  }
  awords.close();
  bwords.close();
  cwords.close();
}

void lpGenerator::writeSigmaFunctions(ofstream* output) {
  int row, triangleNum, identityNum, nonNegativityNum;
  int i, a, b, c, x;

  (*(output)) << "\nSubject To\n";

  row = 0;
  triangleNum = 0;
  identityNum = 0;
  nonNegativityNum = 0;

  /* ******************* */
  /* triangle inequality */
  /* ******************* */
  cout << "starting triangles" << endl << numClasses - 1 << endl << lengthOfWord << endl;
  for (x = 1; x <= numClasses - 1; x++) {
     for (i = 1; i <= lengthOfWord; i++) {
        for (a = 1; a <= ALPHABET_SIZE; a++) {
           for (b = 1; b <= ALPHABET_SIZE; b++) {
              for  (c = a + 1; c <= ALPHABET_SIZE; c++) {
                 if (b != a && b != c) {
                    row++;
                    triangleNum++;
                    /* adheres to the formula -s(X,i,a,b) + -s(X,i,b,c) + s(X,i,a,c) <= 0 */
                    //insertInequality(G, -1, row, mapping->variableToIndex[x][i][a][b]);
                    //insertInequality(G, -1, row, mapping->variableToIndex[x][i][b][c]);
                    //insertInequality(G, 1, row, mapping->variableToIndex[x][i][a][c]);
                    //h[row][1] = 0.0; /* the constant side of the inequality */

                    int order1=residues2[a-1]<residues2[b-1];
                    int order2=residues2[b-1]<residues2[c-1];
                    int order3=residues2[a-1]<residues2[c-1];

                    char Ares = residues2[a-1];
                    char Bres = residues2[b-1];
                    char Cres = residues2[c-1];
                    (*(output)) << " t" << triangleNum << ": - s" << x << "_" << i << "_" << (char) ((order1)?(Ares):(Bres)) <<(char) ((order1)?(Bres):(Ares));
                    (*(output)) << " - s" << x << "_" << i << "_" << (char) ((order2)?(Bres):(Cres)) <<(char) ((order2)?(Cres):(Bres));
                    (*(output)) << " + s" << x << "_" << i << "_" <<( char) ((order3)?(Ares):(Cres)) <<(char) ((order3)?(Cres):(Ares));
                    (*(output)) << " <= 0.0\n";
                 }
              }
           }
        }
     }
  }

  /*TODO: make epsilon = 1 */
  /* ******************* */
  /* identity inequality */
  /* ******************* */
  for (x = 1; x <= numClasses - 1; x++) {
     for (i = 1; i <= lengthOfWord; i++) {
        for (a = 1; a <= ALPHABET_SIZE; a++) {
           for (b = a + 1; b <= ALPHABET_SIZE; b++) {
              row++;
              identityNum++;
              /* adheres to the formula s(X,i,a,a) + -s(X,i,a,b) <= -(epsilon) */
              //insertInequality(G, 1, row, mapping->variableToIndex[x][i][a][a]);
              //insertInequality(G, -1, row, mapping->variableToIndex[x][i][a][b]);
              //h[row][1] = -(EPSILON);

              int order=residues2[a-1]<residues2[b-1];

              char Ares = residues2[a-1];
              char Bres = residues2[b-1];

              (*(output)) << " i" << identityNum << ": s" << x << "_" << i << "_" << residues2[a-1] << residues2[a-1];
              (*(output)) << " - s" << x << "_" << i << "_" <<(char) ((order)?(Ares):(Bres)) <<(char)((order)?(Bres):(Ares));
              (*(output)) << " <= " << EPSILON << endl;

              row++;
              identityNum++;
              //insertInequality(G, 1, row, mapping->variableToIndex[x][i][b][b]);
              //insertInequality(G, -1, row, mapping->variableToIndex[x][i][a][b]);

              //h[row][1] = -(EPSILON);
              (*(output)) << " i" << identityNum << ": s" << x << "_" << i << "_" << residues2[b-1] << residues2[b-1];
              (*(output)) << " - s" << x << "_" << i << "_" << (char)((order)?(Ares):(Bres)) << (char)((order)?(Bres):(Ares));
              (*(output)) << " <= " << EPSILON << endl;
           }
        }
     }
  }

  /* ************************* */
  /* non-negativity inequality */
  /* ************************* */
  //it is very possible that we are missing many of the non-negativity inequalities for the errors in here.
  for (x = 1; x <= numClasses - 1; x++) {
     for (i = 1; i <= lengthOfWord; i++) {
        for (a = 1; a <= ALPHABET_SIZE; a++) {
           for (b = a; b <= ALPHABET_SIZE; b++) {
              row++;
              nonNegativityNum++;
              /* adheres to the formula -s(X,i,a,b) <= 0 */
              //insertInequality(G, -1, row, mapping->variableToIndex[x][i][a][b]);
              //h[row][1] = 0.0; /* the constant side of the inequality */
              int order = residues2[a-1]<residues2[b-1];
              char Ares = residues2[a-1];
              char Bres = residues2[b-1];

              (*(output)) << " n" << nonNegativityNum << ": - s" << x << "_" << i << "_" << (char)((order)?(Ares):(Bres)) << (char)((order)?(Bres):(Ares));
              (*(output)) <<string( " <= 0.0\n");
           }
        }
     }
  }
}

void lpGenerator::writeErrorFunctions(string alphaTSname, string betaTSname, string coilTSname, ofstream* output) {

  /*Process k closest A from the A class*/
  /*d_wv >= distance(w,v) - t_a*/
  int a_querynum = TRSSize;
  int b_querynum = TRSSize;
  int c_querynum = TRSSize;

  string queryWord;
  string nearestNeighbor;
  string alphanearestNeighbor;
  string betanearestNeighbor;

  if (TRSSize > -1) {
     a_querynum = TRSSize;
     b_querynum = TRSSize;
     c_querynum = TRSSize;
  }

  int i;
  int j;
  int m;
  int aa_count = 0;
  ifstream alphaTS(alphaTSname.c_str());
  ifstream betaTS(betaTSname.c_str());
  ifstream coilTS(coilTSname.c_str());
  ifstream alphaTS2(alphaTSname.c_str());
  ifstream betaTS2(betaTSname.c_str());
  ifstream coilTS2(coilTSname.c_str());
  for(i = 0; i < 3; i++) {
	  alphaTS >> queryWord;
	  alphaTS2 >> queryWord;
	  betaTS >> queryWord;
	  betaTS2 >> queryWord;
	  coilTS >> queryWord;
	  coilTS2 >> queryWord;
  }
  for (i = 0; i < TRSSize; i++) {


     alphaTS >> queryWord;
     for (j = 1; j <= k; j++) {
        aa_count++;
        nearestNeighbor = (char*)AANeighbors[i][j];
        (*(output)) << " d_aa_nn_" << aa_count << ": - d_" << queryWord << "_" << j -1 << " <= 0.0\n";
        (*(output)) << " d_aa_" << aa_count << ": - d_" << queryWord << "_" << j -1;

        for (m = 0; m < lengthOfWord; m++) {
           int order = queryWord[m] < nearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = nearestNeighbor[m];
           (*(output)) << " + s" << 1 << "_" << m + 1 << "_" <<(char) ((order)?(Qres):(NNres)) << (char)((order)?(NNres):(Qres));
        }
        (*(output)) << " - t_a <= 0.0\n";
     }
  }

  /*Process l closest B from the A class*/
  int ab_count = 0;
    queryWord = "";
  for (i = 0; i < a_querynum; i++) {
     alphaTS2 >> queryWord;
     for (j = 0; j < l; j++) {
        ab_count++;
        nearestNeighbor =(char*) ABNeighbors[i][j];

        if (j == 0)
          (*(output)) << " d_ab_nn_" << ab_count << ": - d_" << queryWord << " <= 0.0\n";

        (*(output)) << " d_ab_" << ab_count << ": - d_" << queryWord;

        for (m = 0; m < lengthOfWord; m++) {
           int order = queryWord[m] < nearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = nearestNeighbor[m];

           (*(output)) << " - s" << 2 << "_" << m + 1 << "_" << (char)((order)?(Qres):(NNres)) << (char)((order)?(NNres):(Qres));
        }

        (*(output)) << " + t_b <= -1.0\n";
     }
     //this code processes daa < dab
     for (j = 0; j < l; j++) {
        ab_count++;
        alphanearestNeighbor = (char*) AANeighbors[i][1];
        betanearestNeighbor =(char*) ABNeighbors[i][j];

        if (j == 0)
          (*(output)) << " d_ab_nnf_" << ab_count << ": - d_" << queryWord << " <= 0.0\n";

        (*(output)) << " d_abf_" << ab_count << ": - d_" << queryWord;

        for (m = 0; m < lengthOfWord; m++) {
           int order = queryWord[m] < alphanearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = alphanearestNeighbor[m];

           (*(output)) << " + s" << 1 << "_" << m + 1 << "_" << (char)((order)?(Qres):(NNres)) << (char)((order)?(NNres):(Qres));
        }
        for (m = 0; m < lengthOfWord; m++) {
           int order = queryWord[m] < betanearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = betanearestNeighbor[m];

           (*(output)) << " - s" << 2 << "_" << m + 1 << "_" << (char)((order)?(Qres):(NNres)) << (char)((order)?(NNres):(Qres));
        }

        (*(output)) << "  <= 0.0\n";
     }

  }

  /*Process k closest B from the B class*/
  int bb_count = 0;
  for (i = 0; i < b_querynum; i++) {
    betaTS >> queryWord;
     for (j = 1; j <=  k; j++) {
        bb_count++;
        nearestNeighbor =(char*)BBNeighbors[i][j];
        (*(output)) << " d_bb_nn_" << bb_count << ": - d_" << queryWord << "_" << j - 1 << " <= 0.0\n";
        (*(output)) << " d_bb_" << bb_count << ": - d_" << queryWord << "_" << j - 1;

        for (m = 0; m < lengthOfWord; m++) {
           int order=queryWord[m]<nearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = nearestNeighbor[m];

           (*(output)) << " + s" << 2 << "_" << m + 1 << "_" << (char)((order)?(Qres):(NNres)) << (char)((order)?(NNres):(Qres));
        }

        (*(output)) << " - t_b <= 0.0\n";
     }
  }
  printf("%d b-inclass error inequalities\n", bb_count);

  /*Process l closest A from the B class*/
  int ba_count = 0;
  for (i = 0; i < b_querynum; i++) {
    betaTS2 >> queryWord;
     for (j = 0; j < l; j++) {
        ba_count++;
        nearestNeighbor = (char*) BANeighbors[i][j];
        if (j == 0)
        (*(output)) << " d_ba_nn_" << ba_count << ": - d_" << queryWord << " <= 0.0\n";

      (*(output)) << " d_ba_" << ba_count << ": - d_" << queryWord;

        for (m = 0; m < lengthOfWord; m++) {
           int order=queryWord[m]<nearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = nearestNeighbor[m];

           (*(output)) << " - s" << 1 << "_" << m + 1 << "_" << (char)((order)?(Qres):(NNres)) << (char)((order)?(NNres):(Qres));
        }

        (*(output)) << " + t_a   <= -1.0\n";

     }
     //this code processes daa < dab
     for (j = 0; j < l; j++) {
        ab_count++;
        alphanearestNeighbor = (char*) BANeighbors[i][j];
        betanearestNeighbor =(char*) BBNeighbors[i][1];

        if (j == 0)
          (*(output)) << " d_ba_nnf_" << ab_count << ": - d_" << queryWord << " <= 0.0\n";

        (*(output)) << " d_baf_" << ab_count << ": - d_" << queryWord;

        for (m = 0; m < lengthOfWord; m++) {
           int order = queryWord[m] < betanearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = betanearestNeighbor[m];

           (*(output)) << " + s" << 2 << "_" << m + 1 << "_" << (char)((order)?(Qres):(NNres)) << (char)((order)?(NNres):(Qres));
        }
        for (m = 0; m < lengthOfWord; m++) {
           int order = queryWord[m] < alphanearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = alphanearestNeighbor[m];

           (*(output)) << " - s" << 1 << "_" << m + 1 << "_" << (char)((order)?(Qres):(NNres)) << (char)((order)?(NNres):(Qres));
        }

        (*(output)) << "  <= 0.0\n";
     }
  }
  printf("%d b-imposter error inequalities\n", ba_count);

  /*Process l closest A from the C class*/
  int ca_count = 0;
  for (i = 0; i < c_querynum; i++) {
    coilTS >> queryWord;
     for (j = 0; j < l; j++) {
        ca_count++;
        nearestNeighbor =(char*) ACNeighbors[i][j];
        if (j == 0)
          (*(output)) << " d_ca_nn_" << ca_count << ": - d_" << queryWord << "_A <= 0.0\n";

        (*(output)) << " d_ca_" << ca_count << ": - d_" << queryWord << "_A";

        for (m = 0; m < lengthOfWord; m++) {
           int order=queryWord[m]<nearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = nearestNeighbor[m];

           (*(output)) << " - s" << 1 << "_" << m + 1 << "_" << (char)((order)?(Qres):(NNres)) << (char)((order)?(NNres):(Qres));
        }

        (*(output)) << " + t_a   <= -1.0\n";
     }
  }
  printf("%d c-a_imposter error inequalities\n", ca_count);

  /*Process l closest B from the C class*/
  int cb_count = 0;
  for (i = 0; i < c_querynum; i++) {
    coilTS2 >> queryWord;
     for (j = 0; j < l; j++) {
        cb_count++;
        nearestNeighbor =(char*) BCNeighbors[i][j];
        if (j == 0)
          (*(output)) << " d_cb_nn_" << cb_count << ": - d_" << queryWord << "_B <= 0.0\n";

        (*(output)) << " d_cb_" << cb_count << ": - d_" << queryWord << "_B";

        for (m = 0; m < lengthOfWord; m++) {
           int order=queryWord[m]<nearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = nearestNeighbor[m];

           (*(output)) << " - s" << 2 << "_" << m + 1 << "_" << (char) ((order)?(Qres):(NNres)) << (char) ((order)?(NNres):(Qres));
        }

        (*(output)) << " + t_b  <= -1.0\n";

     }
  }
  printf("%d c-b_imposter error inequalities\n", cb_count);

  /*SET BOTH TAU'S EQUAL TO EACH OTHER -- COMMENT OUT IF UNNEEDED*/
  // (*(output)) << " tau1: t_a - t_b <= 0.0\n";
  // (*(output)) << " tau2: t_b - t_a <= 0.0\n";
  alphaTS.close();
  alphaTS2.close();
  betaTS.close();
  betaTS2.close();
  coilTS.close();
  coilTS2.close();
  (*(output)) << "End";
}

void lpGenerator::writeErrorFunctionsFromFiles(ofstream* output) {

  /*Process k closest A from the A class*/
  int a_querynum = TRSSize;
  int b_querynum = TRSSize;
  int c_querynum = TRSSize;
  int lengthOfWord = 23;


  if (TRSSize > -1) {
     a_querynum = TRSSize;
     b_querynum = TRSSize;
     c_querynum = TRSSize;
  }

  int aa_count = 0;
  ifstream aafile("nnsearches/15aaneighbors.txt");
  ifstream bbfile("nnsearches/15bbneighbors.txt");
  ifstream bafile("nnsearches/15baneighbors.txt");
  ifstream abfile("nnsearches/15abneighbors.txt");
  ifstream acfile("nnsearches/15acneighbors.txt");
  ifstream bcfile("nnsearches/15bcneighbors.txt");
  ifstream aa2file;
  ifstream ba2file;
  ifstream bb2file;
  ifstream ab2file;

  int ab_count = 0;

  for (int i = 0; i < a_querynum; i++) {


     string queryWord;
     aafile >> queryWord;
     for (int j = 1; j <= k; j++) {
        aa_count++;
	string nearestNeighbor;
        aafile >> nearestNeighbor;
        (*(output)) << "- d_" << queryWord << "_" << j -1 << " <= 0.0\n";
        (*(output)) <<  "- d_" << queryWord << "_" << j -1;

        for (int m = 0; m < lengthOfWord; m++) {
           int order = queryWord[m] < nearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = nearestNeighbor[m];
	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));
           (*(output)) << " + s" << 1 << "_" << m + 1 << "_" <<firstpart << secondpart;
        }
        (*(output)) << " - t_a <= 0.0\n";
     }

  /*Process l closest B from the A class*/
     for (int j = 0; j < l; j++) {
        ab_count++;
	string nearestNeighbor;
        bafile >> nearestNeighbor;

        if (j == 0)
          (*(output)) << "- d_" << queryWord << " <= 0.0\n";

        (*(output)) << "- d_" << queryWord;

        for (int m = 0; m < lengthOfWord; m++) {
           int order = queryWord[m] < nearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = nearestNeighbor[m];

	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));

           (*(output)) << " - s" << 2 << "_" << m + 1 << "_" << firstpart << secondpart;
        }

        (*(output)) << " + t_b <= -1.0\n";
     }
  }
  bafile.close();
  aafile.close();
  aa2file.open("nnsearches/15aaneighbors.txt");
  ba2file.open("nnsearches/15baneighbors.txt");
  ab_count = 0;
  for (int i = 0; i < a_querynum; i++) {
     string queryWord;
     aa2file >> queryWord;
       string alphanearestNeighbor1;
       aa2file >> alphanearestNeighbor1;
       string alphanearestNeighbor2;
       aa2file >> alphanearestNeighbor2;
       for (int j = 0; j < l; j++) {
          ab_count++;
	  string betanearestNeighbor;
          ba2file >> betanearestNeighbor;

          if (j == 0)
            (*(output)) << "- d_" << queryWord << "_F0"  << " <= 0.0\n";

          (*(output)) << "- d_" << queryWord << "_F0";

          for (int m = 0; m < lengthOfWord; m++) {
             int order = queryWord[m] < betanearestNeighbor[m];
             char Qres = queryWord[m];
             char NNres = betanearestNeighbor[m];
	     
	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));


             (*(output)) << " + s" << 2 << "_" << m + 1 << "_" << firstpart << secondpart;
          }
          for (int m = 0; m < lengthOfWord; m++) {
             int order = queryWord[m] < alphanearestNeighbor1[m];
             char Qres = queryWord[m];
             char NNres = alphanearestNeighbor1[m];
	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));

             (*(output)) << " - s" << 1 << "_" << m + 1 << "_" << firstpart << secondpart;
          }

          (*(output)) << " <= 0.0\n";
          if (j == 0)
            (*(output)) << "- d_" << queryWord << "_F1"  << " <= 0.0\n";

          (*(output)) << "- d_" << queryWord << "_F1";

          for (int m = 0; m < lengthOfWord; m++) {
             int order = queryWord[m] < betanearestNeighbor[m];
             char Qres = queryWord[m];
             char NNres = betanearestNeighbor[m];
	     
	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));


             (*(output)) << " + s" << 2 << "_" << m + 1 << "_" << firstpart << secondpart;
          }
          for (int m = 0; m < lengthOfWord; m++) {
             int order = queryWord[m] < alphanearestNeighbor2[m];
             char Qres = queryWord[m];
             char NNres = alphanearestNeighbor2[m];
	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));

             (*(output)) << " - s" << 1 << "_" << m + 1 << "_" << firstpart << secondpart;
          }

          (*(output)) << " <= 0.0\n";
       
     }
  }

  /*Process k closest B from the B class*/
  int bb_count = 0;
  int ba_count = 0;
  for (int i = 0; i < b_querynum; i++) {
    string queryWord;
    bbfile >> queryWord;
     for (int j = 1; j <=  k; j++) {
        bb_count++;
	string nearestNeighbor;
        bbfile >> nearestNeighbor;
        (*(output)) << "- d_" << queryWord << "_" << j - 1 << " <= 0.0\n";
        (*(output)) << "- d_" << queryWord << "_" << j - 1;

        for (int m = 0; m < lengthOfWord; m++) {
           int order=queryWord[m]<nearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = nearestNeighbor[m];
	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));

           (*(output)) << " + s" << 2 << "_" << m + 1 << "_" << firstpart << secondpart;
        }

        (*(output)) << " - t_b <= 0.0\n";
     }

  /*Process l closest A from the B class*/
     for (int j = 0; j < l; j++) {
        ba_count++;
	string nearestNeighbor;
        abfile >> nearestNeighbor;
        if (j == 0)
        (*(output)) << "- d_" << queryWord << " <= 0.0\n";

      (*(output)) << "- d_" << queryWord;

        for (int m = 0; m < lengthOfWord; m++) {
           int order=queryWord[m]<nearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = nearestNeighbor[m];
	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));

           (*(output)) << " - s" << 1 << "_" << m + 1 << "_" << firstpart << secondpart;
        }

        (*(output)) << " + t_a <= -1.0\n";

     }
  
  }
  abfile.close();
  bbfile.close();
  bb2file.open("nnsearches/15bbneighbors.txt");
  ab2file.open("nnsearches/15abneighbors.txt");
  ba_count = 0;
  for (int i = 0; i < b_querynum; i++) {
     string queryWord;
     bb2file >> queryWord;
       string betanearestNeighbor1;
       bb2file >> betanearestNeighbor1;
       string betanearestNeighbor2;
       bb2file >> betanearestNeighbor2;
       for (int j = 0; j < l; j++) {
          ba_count++;
	  string alphanearestNeighbor;
          ab2file >> alphanearestNeighbor;

          if (j == 0)
            (*(output)) <<  "- d_" << queryWord << "_F0" << " <= 0.0\n";

          (*(output)) <<  "- d_" << queryWord << "_F0" ;

          for (int m = 0; m < lengthOfWord; m++) {
             int order = queryWord[m] < alphanearestNeighbor[m];
             char Qres = queryWord[m];
             char NNres = alphanearestNeighbor[m];
	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));

             (*(output)) << " + s" << 2 << "_" << m + 1 << "_" << firstpart << secondpart;
          }
          for (int m = 0; m < lengthOfWord; m++) {
             int order = queryWord[m] < betanearestNeighbor1[m];
             char Qres = queryWord[m];
             char NNres = betanearestNeighbor1[m];
	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));

             (*(output)) << " - s" << 1 << "_" << m + 1 << "_" << firstpart << secondpart;
          }

          (*(output)) << " <= 0.0\n";
          if (j == 0)
            (*(output)) <<  "- d_" << queryWord << "_F1" << " <= 0.0\n";

          (*(output)) <<  "- d_" << queryWord << "_F1" ;

          for (int m = 0; m < lengthOfWord; m++) {
             int order = queryWord[m] < alphanearestNeighbor[m];
             char Qres = queryWord[m];
             char NNres = alphanearestNeighbor[m];
	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));

             (*(output)) << " + s" << 2 << "_" << m + 1 << "_" << firstpart << secondpart;
          }
          for (int m = 0; m < lengthOfWord; m++) {
             int order = queryWord[m] < betanearestNeighbor2[m];
             char Qres = queryWord[m];
             char NNres = betanearestNeighbor2[m];
	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));

             (*(output)) << " - s" << 1 << "_" << m + 1 << "_" << firstpart << secondpart;
          }

          (*(output)) << " <= 0.0\n";
       
   
     }
  }
  printf("%d b-imposter error inequalities\n", ba_count);

  /*Process l closest A from the C class*/
  int ca_count = 0;
  for (int i = 0; i < c_querynum; i++) {
    string queryWord;
    acfile >> queryWord;
     for (int j = 0; j < l; j++) {
        ca_count++;
	string nearestNeighbor;
        acfile >> nearestNeighbor;
        if (j == 0)
          (*(output)) <<  "- d_" << queryWord << "_A <= 0.0\n";

        (*(output)) <<  "- d_" << queryWord << "_A";

        for (int m = 0; m < lengthOfWord; m++) {
           int order=queryWord[m]<nearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = nearestNeighbor[m];
	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));

           (*(output)) << " - s" << 1 << "_" << m + 1 << "_" << firstpart << secondpart;
        }

        (*(output)) << " + t_a <= -1.0\n";
     }
  }
  printf("%d c-a_imposter error inequalities\n", ca_count);

  /*Process l closest B from the C class*/
  int cb_count = 0;
  for (int i = 0; i < c_querynum; i++) {
    string queryWord; 
    bcfile >> queryWord;
     for (int j = 0; j < l; j++) {
        cb_count++;
	string nearestNeighbor;
        bcfile >> nearestNeighbor;
        if (j == 0)
          (*(output)) <<  "- d_" << queryWord << "_B <= 0.0\n";

        (*(output)) <<  "- d_" << queryWord << "_B";

        for (int m = 0; m < lengthOfWord; m++) {
           int order=queryWord[m]<nearestNeighbor[m];
           char Qres = queryWord[m];
           char NNres = nearestNeighbor[m];
	   char firstpart = ( (order) ? (Qres) : (NNres));
	   char secondpart = ( (order) ? (NNres) : (Qres));

           (*(output)) << " - s" << 2 << "_" << m + 1 << "_" <<firstpart << secondpart;
        }

        (*(output)) << " + t_b <= -1.0\n";

     }
  }
  printf("%d c-b_imposter error inequalities\n", cb_count);

  /*SET BOTH TAU'S EQUAL TO EACH OTHER -- COMMENT OUT IF UNNEEDED*/
  // (*(output)) << " tau1: t_a - t_b <= 0.0\n";
  // (*(output)) << " tau2: t_b - t_a <= 0.0\n";
  (*(output)) << "End";
  aafile.close();
  abfile.close();
  bafile.close();
  bbfile.close();
  acfile.close();
  bcfile.close();

}


void lpGenerator::deleteMapping(){
  Mapping m = *(this->mapping);
  delete[] m.indexToClass;
  delete[] m.indexToPosition;
  delete[] m.indexToFirstLetter;
  delete[] m.indexToSecondLetter;

  for(int i = 1; i <= numClasses; i++) {
    for(int j = 1;j <= lengthOfWord; j++) {
      for(int k = 1; k <= ALPHABET_SIZE; k++) {
        delete[] m.variableToIndex[i][j][k];
      }
    }
  }

  for(int i = 1; i <= numClasses; i++) {
    for(int j = 1;j <= lengthOfWord; j++) {
      delete[] m.variableToIndex[i][j];
    }
  }

  for( int i = 1; i <= numClasses; i++) {
    delete[] m.variableToIndex[i];
  }

  delete[] m.variableToIndex;

}

float getdist(int rand,char* word1, char* word2){
  //string weightsfilename = "weights1.txt" ;
  //readinweights(this->lengthOfWord,weightsfilename.c_str());

    //cout << WeightedWordDist2((Pointer) word1,(Pointer) word2);
    //cout << WeightedWordDist((Pointer) word1,(Pointer) word2);
    if (rand != 0){
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
        setdfs(adf,bdf);
        int cset = setcoildf(cdf);
        return wordDist2Beta((Pointer) word1, (Pointer) word2);
    } else {
        return WeightedWordDist2((Pointer) word1,(Pointer) word2);
    //alphadf = adf;
}
}

void writedistfile(string infname, string outname){
    ifstream inf(infname);
    ofstream outf(outname);
    string line;
    int rand = 23;
    string weightsfilename = "weights/weights" + to_string(rand) + ".txt";
    readinweights(23,weightsfilename.c_str());
    if (inf.is_open()){
        while( getline (inf,line)){
            istringstream ss(line);
            string queryword,target;
            ss >> queryword;
            ss >> target;
            char tword[target.size() + 1];
            strcpy(tword,target.c_str());
            char qword[target.size() + 1];
            strcpy(qword,queryword.c_str());
            float dist = getdist(0,qword,tword);
            outf << queryword << " " <<  target << " " << dist << endl;
        }
    }
}


void lpGenerator::gettargetsdist(int nnsearchnum){
    writedistfile("nnsearches5/aaneighbors.txt","nnsearches5/naaneighbors.txt");
    writedistfile("nnsearches5/bbneighbors.txt","nnsearches5/nbbneighbors.txt");
    writedistfile("nnsearches5/ccneighbors.txt","nnsearches5/nccneighbors.txt");
    writedistfile("nnsearches6/aaneighbors.txt","nnsearches6/naaneighbors.txt");
    writedistfile("nnsearches6/bbneighbors.txt","nnsearches6/nbbneighbors.txt");
    writedistfile("nnsearches6/ccneighbors.txt","nnsearches6/nccneighbors.txt");
}


void lpGenerator::getnewtargets(int rand){
    for (int i =0;i < 100;i++){
        ifstream inf ("predictions/23sstrainingsample/predictknn" + to_string(i) + ".txt");
        ofstream outf("predictions/23sstrainingsample/newpredictknn" + to_string(i) + ".txt");
        string line;
        if (inf.is_open()){
            while( getline (inf,line)){
                //process each line here.
                istringstream ss(line);
                string queryword,target;
                ss >> queryword;
                char qword[queryword.size() + 1];
                strcpy(qword,queryword.c_str());
                outf << queryword << " ";

                ss >> target;
                outf << target << " ";
                char tword[target.size() + 1];
                strcpy(tword,target.c_str());
                float dist = getdist(rand,qword,tword);
                outf << dist << " ";

                ss >> target;
                ss >> target;
                ss >> target;
                outf << target << " ";
                strcpy(tword,target.c_str());
                dist = getdist(rand,qword,tword);
                outf << dist << " ";

                ss >> target;
                ss >> target;
                ss >> target;
                outf << target << " ";
                strcpy(tword,target.c_str());
                dist = getdist(rand,qword,tword);
                outf << dist << "\n";
            }
            outf.close();
            inf.close();
        }
    }




}

void lpGenerator::testDF(int rand,char* word1, char* word2){
  //string weightsfilename = "weights1.txt" ;
  //readinweights(this->lengthOfWord,weightsfilename.c_str());

    //cout << WeightedWordDist2((Pointer) word1,(Pointer) word2);
    //cout << WeightedWordDist((Pointer) word1,(Pointer) word2);
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
    setdfs(adf,bdf);
    int cset = setcoildf(cdf);
    setdfs(adf,bdf);
    cset = setcoildf(cdf);
    Distance alphadistfunc = SecondaryDist;
    cout << alphadistfunc((Pointer) word1, (Pointer) word2);
    cout << "\n";
    cout << wordDist2Beta((Pointer) word1, (Pointer) word2);
    cout << "\n";
    cout << wordDist2Coil((Pointer) word1, (Pointer) word2);
    cout << "\n";
    //alphadf = adf;
}


/*===========================================================================*/
/*===============================[ ]===============================*/
/*===========================================================================*/


