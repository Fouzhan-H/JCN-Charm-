#include "CompJointContourNet.h"

#include "stdio.h"

//#include <limits>
#include <float.h>

CompJointContourNet::CompJointContourNet(RectilinearGridData *dset, int n, double * bases, double* sws, IdType pointsNr)
  : fragsSlabIdx(0, std::vector<short>(n)) 
  , Fields(std::vector<Field*>(n)), dataset(dset){
  N = n; 
  M_geom = dataset -> getGeometry();  
  M_topo = dataset -> CellComponentNr() - 1;  
  for (int i = 0 ; i < N; i++) 
     Fields[i] = new Field ( *(bases+i), *(sws+i), pointsNr); 
}

CompJointContourNet::~CompJointContourNet(){
  for (int i = 0; i < N; i++)
    delete Fields[i]; 
}

void CompJointContourNet::compJCN(IdType & bndSlabNr, IdType & SlabNr) {
     
  double point [MaxDomainDim];  

  double domRange[2*MaxDomainDim]; 
  double rngRange[2*MaxRangeDim];  

  double coord[MaxDomainDim][MaxRangeDim];
  
  IdType pointIds[MaxDomainDim];
  
  double thresholds[MaxNrSlicesPerPtope];
  IdType top1[MaxNrSlicesPerPtope];
  IdType top2[MaxNrSlicesPerPtope];
  IdType *toFragment   = top1;
  IdType *newFragments = top2;

  int toFragmentNum   = 0;
  int newFragmentsNum = 0;

  int slice; 
 
  IdType frags[MaxNrSlicesPerPtope];
  IdType ends[MaxNrSlicesPerPtope];
   
  IdType cellsNum = dataset-> GetCellsNr(); 
  //read data

  fetchData(); 

  double dataBounds [MaxDomainDim];
  dataset -> GetSpatialBounds(dataBounds);

  geometry.reset(new PolytopeGeometry(M_topo, M_geom, N, dataBounds)); 

  // Main Loop:
  //
  // For each cell in the dataset ...
  for (IdType c = 0; c < cellsNum ; c++ ){     

    //Compute the domain bounds for the current cell. 
    
    // Get the domain coordinate for the first point in the cell. 
    dataset -> GetCellComponent (c, 0, point); 
    
    // Initialize the min/max range to the first point cordinates.  
    for (int j = 0; j < M_geom; j++)   
      domRange[2*j] = domRange[2*j+1]= point[j]; 
    
    // Comapre the bounds with the remaining points of the current cell. 
    for (int i = 1; i < M_topo+1; i++){
      dataset -> GetCellComponent (c, i, point); 
      
      for (int j = 0; j < M_geom; j++){
        domRange[2*j] = std::min(domRange[2*j], point[j]); 
        domRange[2*j+1] = std::max(domRange[2*j+1], point[j]);  
      }
    }
    

    // Prepare the geometry processor to receive new active polytopes.
    geometry->ResetForNextCell(domRange);
    
    //Collect the range values for each point in the simplex, 
    //compute the range bounds in each dimension, and 
    //create the corresponding point in the geometry processor. 
    for (int i =0; i < this->N; i++){
        rngRange[2*i] = DBL_MAX; 
        //rngRange[2*i] = std::numeric_limits<double>::max(); 
        rngRange[2*i+1]= -DBL_MAX; 
        //rngRange[2*i+1]= std::numeric_limits<double>::lowest(); 
    }
    
    double fv;
    IdType rid;  
    for (int v = 0; v < M_topo+1; v++){
      // Get the domain coordinate and range index for each point in the simplex  
      rid = dataset -> GetCellComponent(c, v, point);  
      
      //Collect the sample values and compute the sample ranges
      for (int r = 0; r < this->N; r++){
        coord[v][r] = fv = this->Fields[r]->pointScalars[rid];

        if (fv < rngRange[2*r]) rngRange[2*r] = fv; 
        if (fv > rngRange[2*r+1]) rngRange[2*r+1] = fv; 
      }
      //Add the vertex in the geometry processor  
      pointIds[v] = geometry -> AddCoord(point, coord[v]); 
    }
  

    //  Build a simplex in the Geometry Processor 
    toFragmentNum = 0; 
    toFragment[toFragmentNum++] = geometry -> BuildFromSimplex(pointIds, M_topo); 
    
    
    //  For each Field 
    for (int f=0; f < this->N; f++) { 
      //    compute cutting planes 
          
      thresholds[0] = this-> Fields[f]-> SlabRangeValue(rngRange[2*f]);  
      
      for (slice = 1; thresholds[slice-1] < rngRange[2*f+1]; slice++)
        thresholds[slice] = thresholds[slice-1] + (this->Fields[f]->slabWidth);  


      if (slice == 1 && thresholds[0] <= rngRange[2*f]){
         // Trivial case: the only cutting plane lies at 
         // or before thte cell minimum, i.e. the cell is
         // fully contained with this quanitzation level. 
         // In short the cell is not cut. 
      }
      else {
        //    call processor.split on list of fragments
        newFragmentsNum = 0; 
        for (int p = 0; p < toFragmentNum; p++){
            geometry -> Split(toFragment[p], 
                              thresholds, slice, 
                              f, frags, ends); 

            for (int i = 0; i <= slice; i++){
               if (frags[i]) newFragments[newFragmentsNum++] = frags[i]; 
            }

        } // loop on each fragment   
        
        std::swap<IdType*> (toFragment, newFragments);
        std::swap<int> (toFragmentNum, newFragmentsNum);
      } // Non-trivial case 
            
          
    } // loon on each field

     
     geometry->StoreActivePolytopes(toFragment, toFragmentNum);
     
 
   } //Main loop on cells     



  // Merge Fragments ..
  // Form the slabs by merging adjacent equivalent fragments, 
  // record slab ids and slab values.
  // Record edges
  // Record facets on the boundary   


  IdType fragsNr = geometry -> GetNrStoredPtopes();   
  // Fragments slab values 
  fragsSlabIdx.resize(fragsNr, std::vector<short>(this -> N)); 
  double rv; 
  for (IdType i = 0; i < fragsNr; i++){
    for (int j = 0; j < this->N; j++){
       rv = geometry -> GetStoredRangeComponent(i,j); 
       fragsSlabIdx[i][j] = this->Fields[j]->SlabRangeIndex(rv);
    }
  }
  
  std::cout << "Number of frags: " << fragsNr << " " << std::endl;
 
  IdType cntr1 = 0; 
  IdType cntr2 = 0; 
  // Find equivalen fragments    
  // Create and initialize the Union-Find structure 
  this->UF.reset (new std::vector<IdType>(fragsNr, -1)); 
  const IdType * uptr = geometry -> GetCenterFacets(); 
  const IdType * vptr = uptr + 1; 
  IdType facetsNr = geometry -> GetFacetsNr(); 
  IdType rootU, rootV; 
  bool equal; 
  for (IdType c = 0; c< facetsNr; c++, uptr += 2, vptr += 2){     
     
      if (*vptr < 0 ){
         cntr2 ++;
         continue; 
      }

      equal = true; 
      for (int i = 0 ; equal && i < this->N; i++)
         equal = (fragsSlabIdx[*uptr][i] == fragsSlabIdx[*vptr][i]); 

      if (equal){
        // the two fragments are merged
        rootU = this -> Find(*uptr);
        rootV = this -> Find(*vptr);
        if (rootU != rootV)
          this -> UF->operator[](rootV) = rootU;
      }
      else {
        // there is a new edge 
        cntr1++;
      }
  }

  IdType SlabsNr = 0 ; 
  for (IdType i = 0; i < fragsNr; i++){
    if (this -> UF->operator[](i) < 0){
       // set slab id 
       this -> UF->operator[](i)= -(1 + SlabsNr);
       // Store the slab range values into the output structure 
       SlabsNr++; 
    }
  }
 
  bndSlabNr = cntr2; 
  SlabNr = SlabsNr; 
  std::cout << "number of faces : " << facetsNr << " " << cntr1 << " " << std::endl;

}


void CompJointContourNet::extractJCN( std::vector<std::vector <double>> & slabs 
                                    , std::set<std::pair<IdType, IdType>> & edges 
                                    , double * facetCntrs, long long * ids ) {
  // Count the number of slabs, assign new slab ids 
  // and store slab values into the output structure 
  // Set slabs range values into the output structure 
  IdType fragsNr = geometry -> GetNrStoredPtopes();   
  for (IdType i = 0, s_idx = 0; i < fragsNr; i++){
    if (this -> UF->operator[](i) < 0){
       // Store the slab range values into the output structure 
       for (int j = 0; j < this -> N; j++) 
         slabs[s_idx][j] = this -> Fields[j]-> fieldBase +  
                                  (this -> Fields[j]->slabWidth) * fragsSlabIdx[i][j]; 
       s_idx++; 
    }
  }
  

  // Find edges. 
  const IdType * uptr = geometry -> GetCenterFacets(); 
  const IdType * vptr = uptr + 1; 
  IdType facetsNr = geometry -> GetFacetsNr(); 
  IdType vertU, vertV, u, v; 
  double * p = facetCntrs; 
  long long * ids_itr = ids; 
  for (IdType c = 0; c < facetsNr; c++, uptr+=2, vptr+=2){
    if (*vptr < 0 ) {
       // The facet is on boundary 
       geometry -> GetFacetCenter(c, p); // retrive  domain coordinate of the facet center 
       vertU = this -> UF -> operator[]( this -> Find (*uptr)); 
       * ids_itr = (-vertU - 1); 
       ids_itr++; 
       p += M_geom;  
       continue; 
    }
    vertU = this -> UF -> operator[](this -> Find (*uptr)); 
    vertV = this -> UF -> operator[](this -> Find (*vptr)); 
    if (vertU != vertV) {
       u = -vertU - 1; 
       v = -vertV - 1; 
       
       if (v < u) std::swap(u,v); 
       edges.insert (std::pair<IdType, IdType>(u,v)); 
    }    
  }
 
}  


// Find the root of a value stored in the union-find
// structure, performing path compression as we go.
IdType CompJointContourNet::Find(IdType x)
{
  IdType entry = this->UF-> operator[](x);
  if (entry < 0) 
    return x;
  else {
    IdType root = this->Find(entry);
    this->UF->operator[](x) = root;
    return root;
  }
}

void CompJointContourNet::fetchData(){
  for (int i = 0; i < N; i++)
    dataset -> fetchInput(i, &(Fields[i]->pointScalars[0]) );
}
