#include "MergJointContourNet.h"


MergJointContourNet::MergJointContourNet(IdType no1, std::vector<std::vector<double>> * slabs_set1, IdType no2, double * slabs_set2, short rdim)
  : UF (std::vector<IdType>(no1+no2, -1)){
  this -> slabsNr_gr1 = no1; 
  this -> slabs_gr1 = slabs_set1; 
  this -> slabsNr_gr2 = no2; 
  this -> slabs_gr2 = slabs_set2; 
  this -> range_dim = rdim; 
  slabsNr_MJCN = 0;
  new_edges.clear();
}

void MergJointContourNet::AddAdjacentSlabs(long long * slab_ids, IdType n){
  long long * p1 = slab_ids; 
  long long * p2 = slab_ids + 1; 
  IdType root1, root2; 
  bool eq;

  for (IdType i = 0; i < n; i++, p1 += 2, p2+=2){    
    // If Slab values of *p1 and *p2 are equal, merge them in UF 
    // otherwise add an edge between the two nodes 
    eq = AreEqualSlabs(*p1, *p2);
    if (eq){
      root1 = this -> Find(*p1); 
      root2 = this -> Find((*p2) + slabsNr_gr1);
      if (root1 != root2){     
        this -> UF[root2] = root1; //TODO: assign the smaller ID as the root
      }
    }else  
       new_edges.insert(std::pair<IdType, IdType>(*p1, (*p2) + slabsNr_gr1));  
        
  }//end - iteration on number of faces  
}


void MergJointContourNet::ExtractJCNSlabs(IdType & SlabsNr){ 
  slabsNr_MJCN = 0;  // TODO:check 
  for (IdType i = 0; i < UF.size(); i++){
    if (UF[i] < 0){
      UF[i] = -(1 + slabsNr_MJCN);
      slabsNr_MJCN++;
    }
  }//end - iteration on UF
  SlabsNr = slabsNr_MJCN; 
}


void MergJointContourNet::ExtractJCNGr( std::vector<std::vector<double>> & slabs
                                      , std::set<std::pair<IdType, IdType>> & edges_new
                                      , std::set<std::pair<IdType, IdType>> & edges_gr1
                                      , IdType * edges_gr2, int edgesNr_gr2 ){
  double * p; 
  for (IdType i = 0, j = 0; j < slabsNr_MJCN; i++){
    if (UF[i] < 0){
      p = i < slabsNr_gr1 ? (slabs_gr1 ->operator[](i)).data() : slabs_gr2 + (i - slabsNr_gr1) * range_dim; 
      for (int k = 0; k < range_dim; k++)
        slabs[j][k] = *p++; 
      j++; 
    }
  }

  AddUpdEdges(edges_new, this->new_edges);
  AddUpdEdges(edges_new, edges_gr1);
  AddUpdEdges(edges_new, edges_gr2, edgesNr_gr2);
}

void MergJointContourNet::AddUpdEdges(std::set<std::pair<IdType, IdType>> & edges_new, 
                                      std::set<std::pair<IdType, IdType>> & edges_old){
  IdType u, v; 
  for (std::set<std::pair<IdType, IdType>>::iterator it = edges_old.begin(); it != edges_old.end(); ++it){
    u = UF.operator[](Find(std::get<0>(*it)));
    v = UF.operator[](Find(std::get<1>(*it)));
    if (u != v){
      u = -u - 1; 
      v = -v - 1; 
      if (v < u) std::swap(u, v); 
      edges_new.insert (std::pair<IdType, IdType>(u, v)); 
    }  
  }
}

void MergJointContourNet::AddUpdEdges(std::set<std::pair<IdType, IdType>> & edges_new, 
                                      IdType * edges_old, int edgesNr){
  IdType u, v; 
  IdType * ep = edges_old; 
  for (int i = 0; i < edgesNr; i++){
    u = UF.operator[](Find( slabsNr_gr1 + (*ep++) ));
    v = UF.operator[](Find( slabsNr_gr1 + (*ep++) ));
    if (u != v){
      u = -u - 1; 
      v = -v - 1; 
      if (v < u) std::swap(u, v); 
      edges_new.insert (std::pair<IdType, IdType>(u, v)); 
    }  
  }
}

IdType MergJointContourNet::Find(IdType x){
  IdType entry = this -> UF[x]; 
  if (entry < 0)
    return x; 
  else {
    IdType root = this -> Find(entry);
    this -> UF[x] = root; 
    return root; 
  }
}

IdType * MergJointContourNet::UpdIdMap(){
  for (IdType i = 0; i < UF.size(); i++){
    if (UF[i] >= 0)
      UF[i] = UF[Find(i)];
  }
  return UF.data(); 
}


bool MergJointContourNet::AreEqualSlabs(IdType id1, IdType id2){
  // Check the slab values of the two given slab ids
  double * p1 = (slabs_gr1 -> operator[](id1)).data(); 
  double * p2 = slabs_gr2 + id2 * range_dim; 
  for (short i = 0; i < range_dim; i++, ++p1, ++p2){
    if (!FP_isequal(*p1, *p2))
      return false;   
  }
    
  return true;   
}

