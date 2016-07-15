#ifndef BOUNDSLABS_H
#define BOUNDSLABS_H

#include <vector>
#include <array>

#include <ckmulticast.h>

extern int Xdiv;
extern int Ydiv;
extern int Zdiv;
extern CkGroupID mCastGrpId; 

class BSlabsMsg : public CMessage_BSlabsMsg {
public:
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
  BoundSlabs(CkMigrateMessage * msg); 
  void Populate(int * st, int * sz , int M_geom, CkCallback & cb); 
  void SndBSlabs(GetBSlabsMsg * msg);
  void UpdSlabIDs(UpdIDsMsg * msg);  

private:
  bool expectGetBound(); 
  BoundSlabs_SDAG_CODE
  std::vector <std::array<double, 2>> ranges; 
  int M_geom; 
  int iter; 
  int x_itr;
  int y_itr; 
  int z_itr; 
  int max_iter;
  // TODO : data structure
};

#endif

