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
  BndFacesMsg (long long n)
    : CMessage_BndFacesMsg(), no(n){}
  long long * ids; //TODO IDType
  long long no;  
};

class UpdIDsMsg : public CkMcastBaseMsg, public CMessage_UpdIDsMsg {
public:
  UpdIDsMsg(CkCallback  cb)
    : cb(cb){}; 
  long long * ids; //TODO IDType
  CkCallback  cb; //TODO: better way to handle this? 
};

class TargetBndMsg : public CkMcastBaseMsg, public CMessage_TargetBndMsg {
public:
  TargetBndMsg(short d, bool tp)
    : dim(d), topPlain(tp){};
  short dim; 
  bool topPlain;
};

class SndBndMsg : public CkMcastBaseMsg, public CMessage_SndBndMsg {
public:
  SndBndMsg(short d, bool tp, CkCallback cb)
    : dim(d), topPlain(tp), cb(cb){};
  short dim; 
  bool topPlain;
  CkCallback cb; //TODO: better way to handle this? 
};


class BoundSlabs : public CBase_BoundSlabs {
public:
  BoundSlabs(int M_geom); 
  virtual ~BoundSlabs();
  virtual void pup(PUP::er &p);
  BoundSlabs(CkMigrateMessage * msg); 
  void Populate(BSPopulateMsg * msg); 
  void SndSlabFaces(TargetBndMsg * msg);
  void UpdSlabIDs(UpdIDsMsg * msg);  
private:
  BoundSlabs_SDAG_CODE
  void FindAdjSlabs(SndBndMsg * sb_msg, BndFacesMsg * fs_msg); 
  bool ExpectGetBound(); 
  short FindBorder(double * point); 
//TODO  void AdjSlabs(long long facetsNr, BndFacesMsg * set1, BndFacesMsg * set2,  long long * adj_msg); 
  void UpdateSlabIDs(short idx);

  std::vector <std::pair<double, double>> ranges; 
  int M_geom; 
  int iter; 
  int x_itr;
  int y_itr; 
  int z_itr; 
  int max_iter;
  // TODO : data structure
  std::vector<std::vector<long long>> plains; 
  long long * slab_ids;
  long long slab_nr;  
};

#endif

