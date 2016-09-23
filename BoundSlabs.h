#ifndef BOUNDSLABS_H
#define BOUNDSLABS_H

#include "CommonDef.h"
#include <vector>
#include <array>

#include <ckmulticast.h>
#include <pup_stl.h>

extern int Xdiv;
extern int Ydiv;
extern int Zdiv;
extern CkGroupID mCastGrpId; 

class BSPopulateMsg : public CMessage_BSPopulateMsg {
public:
  BSPopulateMsg(CkCallback cb, long long s_nr, long long f_nr) 
    : CMessage_BSPopulateMsg(), cb(cb), slabNr(s_nr), facetNr(f_nr){};
  int * st; 
  int * sz; 
  double * facetCntrs; 
  long long * slabIDs;  
  long long slabNr; 
  long long facetNr; 
  CkCallback cb; 
};

class BndFacesMsg : public CMessage_BndFacesMsg {
public:
  BndFacesMsg (int n)
    : CMessage_BndFacesMsg(), no(n){}
  double * cntrs;  // a 2D array of size no*3 (domain dimension) 
  long long * ids; //TODO IDType
  int no;  
};

/* TODO delete
class AdjFacesMsg {
public:
  AdjFacesMsg(int n) {
    no = n;  
    fst_ids = new long long [n];
    snd_ids = new long long [n];
  };
  AdjFacesMsg(){};

  void pup(PUP::er &p){
    p|no; 
    if (p.isUnpacking()){
      fst_ids = new long long [no];
      snd_ids = new long long [no];
    }
    PUParray(p, fst_ids, no);
    PUParray(p, snd_ids, no);
  };

  ~AdjFacesMsg(){
     delete fst_ids;
     delete snd_ids;
  }

  long long * fst_ids; // std::vector<long long> fst_ids;
  long long * snd_ids; // std::vector<long long> snd_ids;  
  int no; 
};*/

class UpdIDsMsg : public CkMcastBaseMsg, public CMessage_UpdIDsMsg {
public:
  UpdIDsMsg(CkCallback  cb)
    : CkMcastBaseMsg(), CMessage_UpdIDsMsg(), cb(cb){}; 
  long long * ids; //TODO IDType
  CkCallback  cb; //TODO: better way to handle this? 
};

class TargetBndMsg : public CkMcastBaseMsg, public CMessage_TargetBndMsg {
public:
  TargetBndMsg(short d, bool tp)
    : CkMcastBaseMsg(), CMessage_TargetBndMsg(), dim(d), topPlain(tp){};
  short dim; 
  bool topPlain;
};

class SndBndMsg : public CkMcastBaseMsg, public CMessage_SndBndMsg {
public:
  SndBndMsg(short d, bool tp, CkCallback cb)
    : CkMcastBaseMsg(), CMessage_SndBndMsg(), dim(d), topPlain(tp), cb(cb){};
  short dim; 
  bool topPlain;
  CkCallback cb; //TODO: better way to handle this? 
};


class BoundSlabs : public CBase_BoundSlabs {
  BoundSlabs_SDAG_CODE
public:
  BoundSlabs(int M_geom); 
  ~BoundSlabs();
  BoundSlabs(CkMigrateMessage * msg); 
  void Populate(BSPopulateMsg * msg); 
  void SndSlabFaces(TargetBndMsg * msg);
  void UpdSlabIDs(UpdIDsMsg * msg);  
private:
  void FindAdjSlabs(SndBndMsg * sb_msg, BndFacesMsg * fs_msg); 
  bool ExpectGetBound(); 
  short FindBorder(double * point); 
  void AdjSlabs(long long facetsNr, BndFacesMsg * set1, BndFacesMsg * set2,  long long * adj_msg); 
  void UpdateSlabIDs(short idx);
  std::vector <std::array<double, 2>> ranges; 
  int M_geom; 
  int iter; 
  int x_itr;
  int y_itr; 
  int z_itr; 
  int max_iter;
  // TODO : data structure
  std::vector <BndFacesMsg *> plains;
  long long * slab_ids;
  long long slab_nr;  
};

#endif

