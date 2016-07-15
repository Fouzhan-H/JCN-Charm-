#ifndef MAIN_H
#define MAIN_H

#define MAX_RANGE_DIM 10

int XDim;
int YDim;
int ZDim; 
int RngDim; 
int sws [MAX_RANGE_DIM];
int bss [MAX_RANGE_DIM];
int Xdiv;
int Ydiv;
int Zdiv;
CkGroupID mCastGrpId; 
CProxy_BoundSlabs BSArray; 
CProxy_JointContourNet JCNArray; 

class Main : public CBase_Main {
public:
  Main(CkArgMsg *msg);
  Main(CkMigrateMessage *msg);
  void Done(); // TODO: receive the final JCN
private:
  void GetParams (CkArgMsg * msg); 
};

#endif
