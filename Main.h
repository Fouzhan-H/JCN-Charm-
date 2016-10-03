#ifndef MAIN_H
#define MAIN_H

#include <pup_stl.h>

#define MAX_RANGE_DIM 10

int XDim;
int YDim;
int ZDim; 
int Range_Dim; 
double sws [MAX_RANGE_DIM];
double bss [MAX_RANGE_DIM];
char fns [MAX_RANGE_DIM][100];
int Xdiv;
int Ydiv;
int Zdiv;
CkGroupID mCastGrpId; 
CProxy_BoundSlabs BSArray; 
double startTime; 

class Main : public CBase_Main {
public:
  Main(CkArgMsg *msg);
  Main(CkMigrateMessage *msg);
//  void Done(); // TODO: receive the final JCN
private:
  void GetParams (CkArgMsg * msg); 
//  double startTime; 
};

#endif
