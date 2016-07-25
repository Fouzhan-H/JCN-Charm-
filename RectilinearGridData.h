#ifndef RECTILINEARGRIDDATA_H
#define RECTILINEARGRIDDATA_H 

#include "CommonDef.h"

#include <array>
#include <string>
#include <vector>

class RectilinearGridData{
public: 
  RectilinearGridData(int M, int * start, int* size, const int * glb_dim, int N, char fileNames [][100]);
  IdType GetCellsNr(); 

  // returns the dimension of the spatial domain 
  int getGeometry(){return M_geom;};

  int CellComponentNr(); 

  void GetSpatialBounds(double *); 
  
  //return the J-th component (domain) of the cell i. 
  IdType GetCellComponent (IdType i, int j , double *);

  void fetchInput(int i, float*); 
  void fetchCharInput(int i, float*); 
private: 
  // start index and size in domain 
  std::vector<int> st; 
  std::vector<int> sz; 
  const int * glb_dim;

  int M_geom; //dimension of the spatiol domain
  int N;      //dimension of the range data

  //TODO std::string * fnames;   
  char (* fnames) [100];
  void Simplicate3D (short cs); 
  //Cache last cell computation
  IdType cCell;
  int cell_x; 
  int cell_y; 
  int cell_z; 
  std::array<short, 4> cell4; 
};

#endif
