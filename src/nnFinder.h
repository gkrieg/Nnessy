/* ********************************************************
 * Author: Spencer Krieger
 * Date: June 12st, 2016
 *
 * Contains the Dispersion trees and runs nearest neighbor queries on them.
 * ********************************************************/

#ifndef NNFINDER_H_
#define NNFINDER_H_

#include <string>
extern "C" {
  #include "dispersionTree.h"
  #include "metricspaces.h"
  #include "word.h"
}


using namespace std;



class nnFinder {
private:
  DispersionTree* AlphaDT;
  DispersionTree* AlphaDT2;
  DispersionTree* AlphaDT3;
  DispersionTree* BetaDT;
  DispersionTree* BetaDT2;
  DispersionTree* BetaDT3;
  DispersionTree* CoilDT;
  MetricSpace* MSAlpha;
  MetricSpace* MSBeta;

public:
//void printDFs();
  void makeDispersionTree(FILE *queriedWords,  string type, int num, bool firstIteration, int rand);
  void loadDispersionTree(FILE *queriedWords, string filename, string type, int rand);
  void saveDispersionTree(char* filename, string type);
  void destroyDTree(char type);
  void queryAlpha(int numNeighbors,    MetricSpace *MS, Pointer*** AllNeighbors, int treenum  );
  void queryCoil(int numNeighbors,    MetricSpace *MS, Pointer*** AllNeighbors, int treenum  );
  void queryBeta(int numNeighbors  , MetricSpace *MS, Pointer*** AllNeighbors, int treenum  );
  nnFinder() {
    AlphaDT = NULL;
    AlphaDT2 = NULL;
    AlphaDT3 = NULL;
    BetaDT = NULL;
    BetaDT2 = NULL;
    BetaDT3 = NULL;
    CoilDT = NULL;
    MSAlpha = NULL;
    MSBeta = NULL;
  }
  ~nnFinder(){
    if (AlphaDT != NULL)
      destroyDispersionTree(AlphaDT);
    if (AlphaDT2 != NULL)
      destroyDispersionTree(AlphaDT2);
    if (AlphaDT3 != NULL)
      destroyDispersionTree(AlphaDT3);
    if(BetaDT != NULL)
      destroyDispersionTree(BetaDT);
    if(BetaDT2 != NULL)
      destroyDispersionTree(BetaDT2);
    if(BetaDT3 != NULL)
      destroyDispersionTree(BetaDT3);
    if(CoilDT != NULL)
      destroyDispersionTree(CoilDT);
    if(MSAlpha != NULL) {
      DestroyMetricSpace(MSAlpha);
    }
    if(MSBeta != NULL) {
      DestroyMetricSpace(MSBeta);
    }
  }
};

#endif
