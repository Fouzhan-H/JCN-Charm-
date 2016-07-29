#include "BoundSlabs.decl.h"
#include "BoundSlabs.h"

#include <math.h>

#define MYEPSILON 1.0e-7

BoundSlabs :: BoundSlabs(int M_geom) 
  : M_geom (M_geom)
  , ranges(std::vector <std::array<double, 2>>(M_geom))
  , plains (std::vector<BndFacesMsg *>(2*M_geom)){
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

BoundSlabs :: BoundSlabs(CkMigrateMessage * msg){}

BoundSlabs::~BoundSlabs(){
  for (short i = 0; i < plains.size(); i ++){
    if (plains[i] != NULL) delete plains[i]; 
  }
  if (slab_ids != NULL) delete slab_ids;   
}

void BoundSlabs :: Populate(BSPopulateMsg * msg){
  long long facetNr = msg -> facetNr; 
  for (int i = 0 ; i < M_geom; i++){
    ranges[i][0] = msg -> st[i];
    ranges[i][1] = msg -> st[i] + msg -> sz[i] -1; 
  }

  // Iterate over all facets 
  //   * determine which border they sit on
  //   * Count the number of faces on each border
  //   * Create new message for each border
  //   * Intialize each border with the corresponding face centers and ids
     
  std::vector<int> cntr = std::vector<int>(2*M_geom,0);
  short * plainIdx = new short [facetNr];  
  short * b_idx = plainIdx;  
  double * p = msg -> facetCntrs; 
  for (long long i = 0; i < facetNr; i++, p += M_geom, b_idx++){ 
    *b_idx = FindBorder(p);
    if (*b_idx == -1) CkExit();
    cntr[*b_idx] += 1; 
  }

  for (int i = 0; i < 2*M_geom; i++)   // messages for each border
    plains[i] = new (cntr[i]*M_geom, cntr[i]) BndFacesMsg(cntr[i]);
  
  short n; 
  long long m;
  b_idx = plainIdx; 
  std::vector<long long> f_idx = std::vector<long long>(2*M_geom, 0);
  for (long long i = 0; i < facetNr; i++, b_idx++){ // initialize messages 
    n = *b_idx; 
    m = f_idx[n];
    // copy point and id to the correspondign boundary
    *(plains[n] -> ids + m)  = msg -> slabIDs[i]; 
    for (int j = 0; j < M_geom; j++)
      *(plains[n] -> cntrs + m * M_geom + j) =  msg -> facetCntrs [i*M_geom + j];
    f_idx [n] += 1;
  }

  slab_nr = msg -> slabNr; 
  slab_ids = new long long [slab_nr];
  for (long long i = 0; i < slab_nr; i++)
    slab_ids [i] = i; 


  (msg->cb).send();
  delete msg; 
  BSIterate();
}


void BoundSlabs :: SndSlabFaces(TargetBndMsg * msg){
  short b_idx = 2*(msg ->dim-1) + (msg->topPlain ? 1 : 0); 
  BndFacesMsg * fs_msg = plains[b_idx];
  // Update slab IDs
  long long tmp; 
  for (long long i = 0; i < fs_msg -> no ; i++){
    tmp = fs_msg -> ids[i]; 
    fs_msg -> ids[i] = slab_ids[tmp];    
  }
 
  CkSetRefNum(fs_msg, iter);
  switch (msg->dim){
    case 1: thisProxy(thisIndex.x-1, thisIndex.y, thisIndex.z).RecvSlabFaces(fs_msg); break;
    case 2: thisProxy(thisIndex.x, thisIndex.y-1, thisIndex.z).RecvSlabFaces(fs_msg); break;
    case 3: thisProxy(thisIndex.x, thisIndex.y, thisIndex.z-1).RecvSlabFaces(fs_msg); break;
  }
 
  delete msg;  
  plains[b_idx] = NULL;
}

void BoundSlabs::FindAdjSlabs(SndBndMsg * sb_msg, BndFacesMsg * adjfs_msg){
  short b_idx = 2*(sb_msg ->dim-1) + (sb_msg->topPlain ? 1 : 0); 
  BndFacesMsg * fs_msg = plains[b_idx];
  AdjFacesMsg * adj_msg = new (fs_msg->no, adjfs_msg -> no) AdjFacesMsg();
  //TODO: Find adjacent slabs , iterate over fs_msg and adjfs_msg
  
  CkSectionInfo cookie; 
  CkGetSectionInfo(cookie, sb_msg);
  CkMulticastMgr * mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch(); 
  mCastGrp -> contribute(sizeof(*adj_msg), adj_msg , CkReduction::set, cookie, sb_msg->cb); 
  
  delete sb_msg; 
  delete adjfs_msg;
  delete fs_msg;   // TODO :check  correctness
  plains [b_idx] = NULL;  
}


void BoundSlabs :: UpdSlabIDs(UpdIDsMsg * msg){
  CkSectionInfo upd_cookie;
  CkGetSectionInfo(upd_cookie, msg);
  CkMulticastMgr * mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  mCastGrp -> contribute(0, NULL, CkReduction::nop, upd_cookie, msg -> cb);  //TODO is nop a right type?

  // Update slab id record
/*TODO  long long n;
  for (long long i = 0 ; i < slab_nr ; i++){
    n = slab_ids [i];
    slab_ids[i] = msg -> ids [n]; //TODO : can cause out of bound memory access?
  }
  */
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
    if ( fabs( point[i] - ranges[i][0] )  <  MYEPSILON){
      // point lays on the lower plane perpendicular to i-th dimension 
      return 2*i; 
    }
    if ( fabs ( point[i] - ranges[i][1]) < MYEPSILON){
      // point lays on the higher plane perpendicular to i-th dimension 
      return (2*i+1); 
    }
  }
  return -1; 
}
 

#include "BoundSlabs.def.h"
