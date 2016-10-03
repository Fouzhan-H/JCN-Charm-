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
extern double startTime;

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
  virtual ~JointContourNet();
  JointContourNet(CkMigrateMessage * msg);
  virtual void pup(PUP::er &p);
  void Start();    
  void RunMergeStep();  
private:
  JointContourNet_SDAG_CODE
  void ComputeJCN(/* Needs to recive a wrapper over data files or and index*/);
  void MergeJCNs(JCNGrMsg * rcvdGr, CkReductionMsg * AdjSlabs); 
  void UpdateBorders(int set1_Nr, long long * set1_ids, int set2_Nr, long long * set2_ids); 
  void NextMergerStep();

  void datasetIndex(int dim, int n, int idx, int & st_i, int & dim_i);
  void JCNtoFile (std::string fname);
  int itr_idx;
  int red_dim;
  int red_idx;
  bool fin_dim; 
  bool upd_sec1; 
  bool upd_sec2;
  // local copy of a jcn graph includes 
  //   * a list of the center coordination (in the domain) of all slabs, and
  //   * a set of slab id pairs representing adjacent slabs
  std::vector <std::vector<double>> * slabs;        // std::unique_ptr<std::vector <std::vector<double>>> slabs; 
  std::set <std::pair<long long, long long>> * edges; // std::unique_ptr<std::set <std::pair<long long, long long>>> edges; 

};

#endif 
