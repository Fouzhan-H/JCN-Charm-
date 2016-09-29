#include "BoundSlabs.decl.h"
#include "BoundSlabs.h"

#include <math.h>
#include <algorithm>

bool FacetOrder(std::pair<std::vector<double>, long long> u
               , std::pair<std::vector<double>, long long> v){
  for (int i = 0; i < u.first.size(); i++){
    if (FP_isequal(u.first[i], v.first[i]))
      continue;
    else 
      return u.first[i] < v.first[i];
  }
  return true; 
}

//---------------------------------------------------------

BoundSlabs :: BoundSlabs(int M_geom) 
  : M_geom (M_geom)
  , ranges(std::vector <std::pair<double, double>>(M_geom))
  , plains (std::vector<std::vector<long long>>(2*M_geom)), slab_ids(NULL), slab_nr(0){
  //Initializae iter and max_iter, max_iter is equal to the height of the reduction tree
  iter = 1; 
  x_itr = y_itr = z_itr = 0; 
  bool fin_dim = (Ydiv == 1 && Zdiv == 1)? true : false;
  int kk = 1; 
  while (kk < Xdiv ){
    x_itr++; 
    kk = kk <<1;  
    if (fin_dim && (!(thisIndex.x % kk == 0  || thisIndex.x % kk == kk -1) || thisIndex.x == 0 || thisIndex.x == Xdiv -1) )
      break; 
  }
  kk = 1; 
  fin_dim = (Zdiv == 1) ? true: false;
  while (kk < Ydiv ){
    y_itr++; 
    kk = kk <<1;  
    if ( fin_dim &&  (!(thisIndex.y % kk == 0  || thisIndex.y % kk == kk -1) || thisIndex.y == 0 || thisIndex.y == Ydiv -1) )
      break; 
  }
  kk = 1; 
  fin_dim = true; 
  while (kk < Zdiv ){
    z_itr++; 
    kk = kk <<1;  
    if ( !(thisIndex.z % kk == 0  || thisIndex.z % kk == kk -1) || thisIndex.z == 0 || thisIndex.z == Zdiv -1)
      break; 
  }
  max_iter = x_itr + y_itr + z_itr; 
   
}

BoundSlabs :: BoundSlabs(CkMigrateMessage * msg){
  slab_ids =  NULL;
}

void BoundSlabs::pup(PUP::er &p){
  CBase_BoundSlabs::pup(p);
  __sdag_pup(p); 
  
  p|M_geom; 
  p|iter; 
  p|x_itr;
  p|y_itr;
  p|z_itr;
  p|max_iter;
  p|slab_nr;
  p|plains; 
  p| ranges;
  int has_slabids = (slab_ids != NULL);
  p| has_slabids; 
  if (has_slabids){
    if (p.isUnpacking()) slab_ids = new long long [slab_nr];  
    PUParray(p, slab_ids, slab_nr); 
  }else 
     slab_ids = NULL;
} 

BoundSlabs::~BoundSlabs(){
  if (slab_ids != NULL) delete slab_ids;   
}

void BoundSlabs :: Populate(BSPopulateMsg * msg){
  long long facetNr = msg -> facetNr; 
  for (int i = 0; i < M_geom; i++){
    ranges[i].first = msg -> st[i];
    ranges[i].second = msg -> st[i] + msg -> sz[i] -1; 
  }
 
  std::vector< std::vector<std::pair<std::vector<double>, long long>>> bnd_faces (2*M_geom);
  // Iterate over all facets 
  //   * determine which border they sit on
  //   * Insert the center point and slab id in the correct boundary plain 
  //   * Sort faces on each border baced on their center points 
  //   * Record slab ids of sorted faces in each border 
  short plainIdx;
  std::vector<double> cntr_point (M_geom);
  double * p = msg -> facetCntrs; 
  for (long long i = 0; i < facetNr; i++, p += M_geom){ 
    plainIdx = FindBorder(p);
    if (plainIdx == -1) {
      CkPrintf("BArray (%i, %i, %i),  Nr %d of facetNr: %d Plain was not found %f %f %f slabID: %d \n", thisIndex.x, thisIndex.y, thisIndex.z, i, facetNr,  *p, *(p+1), *(p+2), msg -> slabIDs[i]);
      CkExit();
    }
    for (short j =0; j < M_geom; j++) 
       cntr_point[j] = p[j]; 
    bnd_faces[plainIdx].push_back(std::pair<std::vector<double>, long long> (cntr_point, msg->slabIDs[i]));  
  }

  for (short i = 0; i < 2*M_geom; i++)  
    std::sort(bnd_faces[i].begin(), bnd_faces[i].end() , FacetOrder);

  for (short i = 0; i < 2*M_geom; i++){
     plains[i].resize(bnd_faces[i].size()); 
     for (long long  j =0 ; j < bnd_faces[i].size(); j++)
       plains[i][j] = bnd_faces[i][j].second;
  }
  
  slab_nr = msg -> slabNr; 
  slab_ids = new long long [slab_nr];
  for (long long i = 0; i < slab_nr; i++)
    slab_ids [i] = i; 

  (msg->cb).send();
  delete msg; 
  
  BSIterate();
}

void BoundSlabs :: UpdateSlabIDs (short idx){
  // Update slab IDs
  long long tmp; 
  for (long long i = 0; i < plains[idx].size() ; i++){
    tmp = plains[idx][i]; 
    plains[idx][i] = slab_ids[tmp];    
  }
}

void BoundSlabs :: SndSlabFaces(TargetBndMsg * msg){
  short b_idx = 2*(msg ->dim-1) + (msg->topPlain ? 1 : 0); 

  UpdateSlabIDs(b_idx);                                  // Update slab IDs
  long long no = plains[b_idx].size();
  BndFacesMsg * ids_msg = new (no) BndFacesMsg(no);       // Create and initailaize out message
  for (long long i = 0; i < plains[b_idx].size(); i++)
    *(ids_msg -> ids + i) = plains[b_idx][i];

 
  CkSetRefNum(ids_msg, iter);
  switch (msg->dim){
    case 1: thisProxy(thisIndex.x-1, thisIndex.y, thisIndex.z).RecvSlabFaces(ids_msg); break;
    case 2: thisProxy(thisIndex.x, thisIndex.y-1, thisIndex.z).RecvSlabFaces(ids_msg); break;
    case 3: thisProxy(thisIndex.x, thisIndex.y, thisIndex.z-1).RecvSlabFaces(ids_msg); break;
  }

 
  delete msg;  
  plains[b_idx].resize(0);
}

void BoundSlabs::FindAdjSlabs(SndBndMsg * sb_msg, BndFacesMsg * adjfs_msg){
  short b_idx = 2*(sb_msg ->dim-1) + (sb_msg->topPlain ? 1 : 0);
  // Update ids in fs_msg
  UpdateSlabIDs(b_idx); 
  
  if (plains[b_idx].size() != adjfs_msg -> no)
    printf("Serious Problems\n");

  int msg_sz = 2 * (adjfs_msg -> no) + 1; 
  long long * adj_slabs_msg = new long long [msg_sz];
  adj_slabs_msg[0] = adjfs_msg -> no; 
  //Find adjacent slabs , iterate over fs_msg and adjfis_msg  
  long long * adj_msg_it = adj_slabs_msg + 1; 
  for (long long i = 0; i < plains[b_idx].size(); i++){
    *adj_msg_it++ = plains[b_idx][i];
    *adj_msg_it++ = *(adjfs_msg -> ids + i);
  } 

  
  CkSectionInfo cookie; 
  CkGetSectionInfo(cookie, sb_msg);
  CkMulticastMgr * mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch(); 
  mCastGrp -> contribute(msg_sz * sizeof(long long), adj_slabs_msg, CkReduction::set, cookie, sb_msg->cb); 
  
  delete sb_msg; 
  delete adjfs_msg;  
  plains[b_idx].resize(0); 
  delete adj_slabs_msg;
}


void BoundSlabs :: UpdSlabIDs(UpdIDsMsg * msg){

  CkSectionInfo upd_cookie;
  CkGetSectionInfo(upd_cookie, msg);
  CkMulticastMgr * mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  mCastGrp -> contribute(0, NULL, CkReduction::nop, upd_cookie, msg -> cb);  //TODO is nop a right type?

  // Update slab id record
  long long n, m;
  for (long long i = 0 ; i < slab_nr ; i++){
    n = slab_ids [i];
    m = msg -> ids [n]; 
    slab_ids[i] = -m - 1; 
  } 
  
  delete msg; 
} 

bool BoundSlabs :: ExpectGetBound(){
  int k, half_k, tmp; 
  
  if (iter <= x_itr){
    k = 1 << iter; 
    tmp = thisIndex.x % k; 
  }else if (iter <= x_itr + y_itr) {
    k = 1 << (iter - x_itr);
    tmp = thisIndex.y % k; 
  } else {
    k = 1 << (iter - x_itr - y_itr);
    tmp = thisIndex.z % k; 
  }
  
  half_k = k >> 1; 
  if (tmp == half_k || tmp == half_k -1){ 
     return true;
  }
  return false;  
}

short BoundSlabs::FindBorder(double * point){
  for (int i = 0; i < M_geom; i++){
    if ( fabs( point[i] - ranges[i].first)  <   1.0e-7 ){
      // point lays on the lower plane perpendicular to i-th dimension 
      return 2*i; 
    }
    if ( fabs ( point[i] - ranges[i].second) < 1.0e-7 /*FP_EPSILON*/){
      // point lays on the higher plane perpendicular to i-th dimension 
      return (2*i+1); 
    }
  }
  return -1; 
}

#include "BoundSlabs.def.h"
