#ifndef COMPJOINTCOUNTOURNET_H
#define COMPJOINTCOUNTOURNET_H
#include "RectilinearGridData.h"
#include "PolytopeGeometry.h"

#include <math.h>

#include <string>
#include <vector>
#include <set>
#include <memory>


#define EPSILON 1.0e-10

struct Field {
  Field (double bs, double sw, int sz)
  : fieldBase(bs), slabWidth(sw)
  , pointScalars (std::vector<float>(sz)){};
  double slabWidth; 
  double fieldBase; 
  std::vector<float> pointScalars; 
  
  // Utility function for converting a value v in field f into
  // the corresponding slab coordinate given the quantization
  // settings.
  double SlabRangeValue(double v){
     double div = (v+EPSILON - fieldBase) / slabWidth; 
     return fieldBase + slabWidth*floor(div); 
  }
  
  short SlabRangeIndex(double v){
     double div = (v+EPSILON - fieldBase) / slabWidth; 
     return floor(div); 
  }
}; //Field

class CompJointContourNet {
public:
  
  CompJointContourNet(RectilinearGridData *dset, int n, double * bases, double* sws, IdType pointsNr);

  void compJCN(IdType & bndSlabNr, IdType & SlabNr);
  void extractJCN (std::vector <std::vector<double>> & slabs, std::set<std::pair<IdType, IdType>> & edges, double * facetCntrs, long long * ids); 
  ~CompJointContourNet();
private:
  // Fields   
  std::vector<Field*> Fields; 
  void fetchData (); 

  // instance of dataset (startIndex and  size in domain)
  int M_geom;
  int M_topo;
  int N; 
  RectilinearGridData * dataset; 

  std::unique_ptr< PolytopeGeometry > geometry; 
 
  std::unique_ptr<std::vector<IdType>> UF; 

  std::vector< std::vector< short >> fragsSlabIdx; 

  IdType Find(IdType x); 

};

#endif
