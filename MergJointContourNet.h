#ifndef MERGJOINTCONTOURNET_H
#define MERGJOINTCONTOURNET_H

#include "CommonDef.h"

#include <utility>
#include <vector>
#include <set>

class MergJointContourNet {
public:
  MergJointContourNet(IdType no1, std::vector<std::vector<double>> * slabs_set1, IdType no2, double * slabs_set2, short rdim);
  void AddAdjacentSlabs(long long * slab_ids, IdType n); 
  void ExtractJCNSlabs (IdType & SlabsNr);
  void ExtractJCNGr( std::vector<std::vector<double>> & slabs
                   , std::set<std::pair<IdType, IdType>> & edges_new
                   , std::set<std::pair<IdType, IdType>> & edges_gr1
                   , IdType * edges_gr2, int edgesNr_gr2 );
  IdType * UpdIdMap();  
private:
  bool AreEqualSlabs(IdType, IdType);  
  IdType Find (IdType x);
  void AddUpdEdges(std::set<std::pair<IdType, IdType>> & edges_new, 
                   std::set<std::pair<IdType, IdType>> & Edges_old); 
  void AddUpdEdges(std::set<std::pair<IdType, IdType>> & edges_new, 
                   IdType * Edges_old, int edgesNr); 

  std::vector<IdType> UF; 
  std::set<std::pair<IdType, IdType>> new_edges;
 
  IdType slabsNr_gr1;
  std::vector<std::vector<double>> * slabs_gr1; 
  IdType slabsNr_gr2;
  double * slabs_gr2; 
  IdType slabsNr_MJCN;  
  short range_dim; 
};

#endif
