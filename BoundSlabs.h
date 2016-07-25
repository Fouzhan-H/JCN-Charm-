#ifndef BOUNDSLABS_H
#define BOUNDSLABS_H

#include <vector>
#include <array>

#include <ckmulticast.h>

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


class BSlabsMsg : public CMessage_BSlabsMsg {
public:
  BSlabsMsg (int n)
    : CMessage_BSlabsMsg(), no(n){}
  double * cntrs; // a 2D array of size no*3 (domain dimension) 
  long long * ids; //TODO IDType
  int no;  
};

class UpdIDsMsg : public CkMcastBaseMsg, public CMessage_UpdIDsMsg {
public:
  UpdIDsMsg(CkCallback  cb)
    : CkMcastBaseMsg(), CMessage_UpdIDsMsg(), cb(cb){} 
  int * ids; //TODO IDType
  CkCallback  cb; //TODO: better way to handle this? 
};


class GetBSlabsMsg : public CkMcastBaseMsg, public CMessage_GetBSlabsMsg {
public:
  GetBSlabsMsg(short d, bool tp, CkCallback cb)
    : CkMcastBaseMsg(), CMessage_GetBSlabsMsg(), dim(d), topPlain(tp), cb(cb){};
  short dim; 
  bool topPlain;
  CkCallback cb; //TODO: better way to handle this? 
};

class BoundSlabs : public CBase_BoundSlabs {
public:
  BoundSlabs(int M_geom); 
  ~BoundSlabs();
  BoundSlabs(CkMigrateMessage * msg); 
  void Populate(BSPopulateMsg * msg); 
  void SndBSlabs(GetBSlabsMsg * msg);
  void UpdSlabIDs(UpdIDsMsg * msg);  
private:
  bool expectGetBound(); 
  short FindBorder(double * point); 
  BoundSlabs_SDAG_CODE
  std::vector <std::array<double, 2>> ranges; 
  int M_geom; 
  int iter; 
  int x_itr;
  int y_itr; 
  int z_itr; 
  int max_iter;
  // TODO : data structure
  std::vector <BSlabsMsg *> plains ;
  long long * slab_ids;
  long long slab_nr;  
/*  std::vector <std::vector <std::vector<double>> * > lPlainsFaces; 
  std::vector <std::vector <long long> * > lPlainsIDs; 
  std::vector <std::vector <std::vector<double>> * > uPlainsFaces; 
  std::vector <std::vector <long long> * > uPlainsIDs; */
};

#endif

