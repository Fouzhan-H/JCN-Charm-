#include "BoundSlabs.decl.h"
#include "BoundSlabs.h"


BoundSlabs :: BoundSlabs(int M_geom) 
  : M_geom (M_geom)
  , ranges(std::vector <std::array<double, 2>>(M_geom)){
  
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


void BoundSlabs :: Populate(int * st , int * sz, int M_geom, CkCallback & cb){

  for (int i = 0 ; i < M_geom; i++){
    ranges[i][0] = st[i];
    ranges[i][1] = st[i] + sz[i] -1; 
  }
  
  //TODO: initailaize with Values  
  
  cb.send();
  BSIterate();
}


void BoundSlabs :: SndBSlabs(GetBSlabsMsg * msg){
  CkSectionInfo cookie; 
  CkGetSectionInfo(cookie, msg);
  BSlabsMsg * bs_msg = new (5, 7) BSlabsMsg(); //TODO: intialiaze the message properly 
  CkMulticastMgr * mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  mCastGrp -> contribute(sizeof(*bs_msg), bs_msg , CkReduction::set, cookie, msg->cb); 

  delete msg;  
}

void BoundSlabs :: UpdSlabIDs(UpdIDsMsg * msg){
  CkSectionInfo upd_cookie;
  CkGetSectionInfo(upd_cookie, msg);
  CkMulticastMgr * mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  mCastGrp -> contribute(0, NULL, CkReduction::nop, upd_cookie, msg -> cb);  //TODO is nop a right type?
  // TODO: update the slab ids
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

#include "BoundSlabs.def.h"
