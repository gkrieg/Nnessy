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
#include <iostream>
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

void lpGenerator::buildSaveTrees(char* alpha, char* beta, char* alphaSave, char* betaSave) {
  const char* AlphaInputFile = (const char* ) alpha;
  const char* BetaInputFile = (const char* ) beta;
  string r = "r";
  FILE* alphaFile = fopen(AlphaInputFile, r.c_str());
  FILE* betaFile = fopen(BetaInputFile, r.c_str());
  int rand = 0;
  thread alphatree(&lpGenerator::buildalphatree, this, alphaFile, 1, true, rand);
  thread betatree(&lpGenerator::buildbetatree, this, betaFile, 1, true, rand);
  alphatree.join();
  betatree.join();
  nnfinder1.saveDispersionTree(alphaSave, "alpha");
  nnfinder4.saveDispersionTree(betaSave, "beta");
}
void lpGenerator::printDFs(){
    //nnfinder1.printDFs();
}




void lpGenerator::findNeighborsOneClass(char* TreeIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand) {
  const char* TreeFilename = (const char* ) TreeIF;
  const char* QueryFilename = (const char* ) QueryIF;
  string r = "r";
  FILE* TreeFile = fopen(TreeFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());

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
    if(numneighbors > 10){
        for (int j = 0;j < numneighbors; j++) {
          outfile << (char*) this->AANeighbors[i][j];
          if (j != numneighbors - 1)
              outfile << " ";
        }
    }
    else{

        for (int j = 1;j < numneighbors; j++) {
          outfile << (char*) this->AANeighbors[i][j];
          if (j != numneighbors - 1)
              outfile << " ";
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

void lpGenerator::findNeighborsOneClassCombined(char* TreeIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand, float epsilon, float rho) {
  const char* TreeFilename = (const char* ) TreeIF;
  const char* QueryFilename = (const char* ) QueryIF;
  setEpsilonRho(epsilon,rho);
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
  char* outword = new char[lengthOfWord + 1];
  for (int i = 0;i < TRSSize;i++) {
    infile >> inword;
    strncpy(outword,inword.c_str(),lengthOfWord);
    outfile << outword << " ";
    if(numneighbors > 10){
        for (int j = 0;j < numneighbors; j++) {
            char* combinedword = (char*) this->AANeighbors[i][j];
            strncpy(outword,combinedword,lengthOfWord);
            outword[lengthOfWord] = '\0';
          outfile << outword;
          if (j != numneighbors - 1)
              outfile << " ";
        }
    }
    else{

        for (int j = 1;j < numneighbors; j++) {
            char* combinedword = (char*) this->AANeighbors[i][j];
            strncpy(outword,combinedword,lengthOfWord);
            outword[lengthOfWord] = '\0';
          outfile << outword;
          if (j != numneighbors - 1)
              outfile << " ";
        }
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

void lpGenerator::findDistances(string QueryFilename,string outputfile)
{
    
  ifstream QueryFile(QueryFilename.c_str());
  string** inputwords;
  string line;
  int wordListSize = 0;
  string firstline;
  getline(QueryFile,firstline);
  cout << firstline << endl;
  //firstline.erase(firstline.find("\n"));
  int sizeofline = 0;
  while(getline(QueryFile,line)){
      sizeofline = 1;
      int next = line.find(" ");
      int last = 0;
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

void lpGenerator::findNeighborsOneClassTwoDistance(char* AlphaIF, char* BetaIF, char* QueryIF, int numneighbors, string outputfile, bool firstIteration, int type, int rand) {
  const char* AlphaFilename = (const char* ) AlphaIF;
  const char* BetaFilename = (const char* ) BetaIF;
  const char* QueryFilename = (const char* ) QueryIF;
  string r = "r";
  FILE* AlphaFile = fopen(AlphaFilename, r.c_str());
  FILE* BetaFile = fopen(BetaFilename, r.c_str());
  FILE* QueryFile = fopen(QueryFilename, r.c_str());

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
  string a = outputfile + "a";
  ofstream outfilea(a.c_str());
  for (int i = 0;i < TRSSize;i++) {
      for (int j = 0; j < numneighbors - 1; j++){
          if (type == 1){
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->AANeighbors[i][j + 1] << " " << wordDist2Alpha(MSQuery->points[i], this->AANeighbors[i][j + 1]) << " " ;
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->BBNeighbors[i][j] << " " << wordDist2Beta(MSQuery->points[i], this->BBNeighbors[i][j]) << endl;
            }
          else if (type == 2) {
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->AANeighbors[i][j] << " " << wordDist2Alpha(MSQuery->points[i], this->AANeighbors[i][j]) << " " ;
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->BBNeighbors[i][j+1] << " " << wordDist2Beta(MSQuery->points[i], this->BBNeighbors[i][j+1]) << endl ;
          }
          else if (type == 3) {
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->AANeighbors[i][j] << " " << wordDist2Alpha(MSQuery->points[i], this->AANeighbors[i][j]) << " " ;
          outfilea << (char*) MSQuery->points[i] << " ";
              outfilea <<(char*) this->BBNeighbors[i][j] << " " << wordDist2Beta(MSQuery->points[i], this->BBNeighbors[i][j]) << endl;
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
void lpGenerator::predictProtein(char* AlphaIF, char* BetaIF, char* ProteinIF, double alphatau, double betatau, string filename, int numneighbors, int rand) {

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
  char predictedclass;
  ofstream outfile(filename.c_str());
  for (int i = 0; i < TRSSize; i++) {
	  double alphadist = 0;
      double betadist = 0;
      for (int j = 0;j < numneighbors; j++){
		   alphadist += wordDist2Alpha(MSAlpha->points[i], AANeighbors[i][j]);
		   betadist += wordDist2Beta(MSAlpha->points[i], BBNeighbors[i][j]);
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

/*===========================================================================*/
/*===============================[ ]===============================*/
/*===========================================================================*/

