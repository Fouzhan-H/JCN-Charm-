#include "BoundSlabs.decl.h"
#include "BoundSlabs.h"

#include <math.h>

#define MYEPSILON 1.0e-7

BoundSlabs :: BoundSlabs(int M_geom) 
  : M_geom (M_geom)
  , ranges(std::vector <std::array<double, 2>>(M_geom))
  , plains (std::vector<BSlabsMsg *>(2*M_geom)){
  //TODO: initializae iter and max_iter
  iter = 0; 
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
    plains[i] = new (cntr[i]*M_geom, cntr[i]) BSlabsMsg(cntr[i]);
  
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


void BoundSlabs :: SndBSlabs(GetBSlabsMsg * msg){
  CkSectionInfo cookie; 
  CkGetSectionInfo(cookie, msg);
  short b_idx = 2*(msg ->dim-1) + (msg->topPlain ? 1 : 0); 
  BSlabsMsg * bs_msg = plains[b_idx];
  // Update slab IDs
  long long tmp; 
  for (long long i = 0; i < bs_msg -> no ; i++){
    tmp = bs_msg -> ids[i]; 
    bs_msg -> ids[i] = slab_ids[tmp];    
  }
  CkMulticastMgr * mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  mCastGrp -> contribute(sizeof(*bs_msg), bs_msg , CkReduction::set, cookie, msg->cb); 

  delete msg;  
  plains[b_idx] = NULL;
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

bool BoundSlabs :: expectGetBound(){
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
  if (tmp == half_k || tmp == half_k -1)
    return true;
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
