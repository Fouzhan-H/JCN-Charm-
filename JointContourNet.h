#ifndef JOINTCONTOURNET_H
#define JOINTCONTOURNET_H

#include "BoundSlabs.h"

extern int XDim;
extern int YDim;
extern int ZDim;
extern int Xdiv;
extern int Ydiv;
extern int Zdiv;
extern char fns [][100]; 
extern int Range_Dim; 
extern double sws [];
extern double bss [];
extern CProxy_BoundSlabs BSArray; 
extern CkGroupID mCastGrpId;


class JCNGrMsg : public CMessage_JCNGrMsg{
public:
  double * vexs;  // 2D array of size vexCnt*N (range dimension)
  long long * edgs; //TODO IDType // 2D array of size 2*edgCnt
//  std::vector<std::vector<double>> slabval; 
  int slabCnt; //TODO IDType 
  int edgCnt; //TODO IDtype
};


class JointContourNet : public CBase_JointContourNet{
public:
  JointContourNet();
  ~JointContourNet();
  JointContourNet(CkMigrateMessage * msg);
  void Start(CkCallback & cb);    
  void RunMergeStep();  

private:
  JointContourNet_SDAG_CODE
  void ComputeJCN(/* Needs to recive a wrapper over data files or and index*/);
  void MergeJCNs(JCNGrMsg * sndGr, CkReductionMsg * fstBr); 
  void UnfoldRecvBSlabs(CkReductionMsg * msg); 
  void UpdateBorders(); 
  void NextMergerStep();

  void datasetIndex(int dim, int n, int idx, int & st_i, int & dim_i);
  int itr_idx;
  int red_dim;
  int red_idx;
  bool fin_dim; 
  CkCallback * finish_cb; 
  // local copy of a jcn graph includes 
  //   * a list of the center coordination (in the domain) of all slabs, and
  //   * a set of slab id pairs representing adjacent slabs
  std::vector <std::vector<double>> slabs; 
  std::set <std::pair<long long, long long>> edges; 

};

#endif 
