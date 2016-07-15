#include "JointContourNet.decl.h"
#include "JointContourNet.h"


JointContourNet :: JointContourNet(){
  fin_dim = (Ydiv == 1 && Zdiv == 1) ? true: false; 
  itr_idx = 0;  
  red_dim =1 ; red_idx = 0; 
  NextMergerStep(); 
}
 
JointContourNet :: ~JointContourNet(){
  if (finish_cb != NULL)
    delete finish_cb;
}

void JointContourNet::Start(CkCallback & cb){
  // Compute the corressponding data index 
  finish_cb = new CkCallback (cb); 
  int st_x, dim_x, st_y, dim_y, st_z, dim_z; 

  datasetIndex(XDim, Xdiv, thisIndex.x, st_x, dim_x); 
  datasetIndex(YDim, Ydiv, thisIndex.y, st_y, dim_y); 
  datasetIndex(ZDim, Zdiv, thisIndex.z, st_z, dim_z);

  int st[3] = {st_x, st_y, st_z};
  int sz[3] = {dim_x, dim_y, dim_z}; 

  //TODO: create RecGrid Pointer 
  //TODO: call ComputeJCN 

  ComputeJCN(); 

  //TODO Populate the coressponding BoundSlab Array element   
  CkCallback mrg_cb(CkIndex_JointContourNet::RunMergeStep()
                   , CkArrayIndex3D(thisIndex.x, thisIndex.y, thisIndex.z), thisProxy);
//TODO  BoundSlabs * bslocal = BSArray(thisIndex.x, thisIndex.y, thisIndex.z).ckLocal();
//  if (bslocal == NULL)
     BSArray(thisIndex.x, thisIndex.y, thisIndex.z).Populate(st, sz, 3, mrg_cb);
//  else  
//    bslocal->Populate( st,  sz , 3, mrg_cb);

}  

JointContourNet :: JointContourNet(CkMigrateMessage * msg){}

void JointContourNet::datasetIndex(int dim, int n, int idx, int & st_i, int & dim_i){
  int quot = (dim-1) / n; 
  int rem = (dim-1) % n; 

  if (idx < rem){
    st_i = (quot + 1) * idx;
    dim_i = quot + 2;
  } else {
    st_i = rem * (quot + 1) + (idx - rem)* quot;
    dim_i = quot+1;
  }
}


void JointContourNet :: RunMergeStep(){

  int k = 1 << red_idx; 

  int lvl_id;
  switch (red_dim){
    case 1: lvl_id = thisIndex.x % k;  break; 
    case 2: lvl_id = thisIndex.y % k;  break; 
    case 3: lvl_id = thisIndex.z % k;  break; 
  }

  if (lvl_id == 0){ //Chare Index is a merger
    // call BoundSlabs.GetBoundSlabs with callback to RecvBSlabSet1 and RecvBSlabsSet2
  
    // Add array indices
    CProxySection_BoundSlabs lplain_proxy, uplain_proxy;  
    int tmp;
    switch (red_dim){ // TODO: currently works only for power of 2
      case 1: tmp = thisIndex.x + (k >> 1); 
              lplain_proxy = CProxySection_BoundSlabs::ckNew( BSArray, tmp, tmp, 1, 
                                                              thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1);
              uplain_proxy = CProxySection_BoundSlabs::ckNew( BSArray, tmp-1, tmp-1, 1, 
                                                              thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1);
              break; 
      case 2: tmp = thisIndex.y + (k >> 1); 
              lplain_proxy = CProxySection_BoundSlabs::ckNew( BSArray, 0, Xdiv-1, 1, 
                                                              tmp, tmp, 1, thisIndex.z, thisIndex.z, 1);
              uplain_proxy = CProxySection_BoundSlabs::ckNew( BSArray, 0, Xdiv-1, 1, 
                                                              tmp-1, tmp-1, 1, thisIndex.z, thisIndex.z, 1);
              break; 
      case 3: tmp = thisIndex.z + (k >> 1); 
              lplain_proxy = CProxySection_BoundSlabs::ckNew( BSArray, 0, Xdiv-1, 1, 
                                                                       0, Ydiv-1, 1, tmp, tmp, 1);
              uplain_proxy = CProxySection_BoundSlabs::ckNew( BSArray, 0, Xdiv-1, 1, 
                                                                       0, Ydiv-1, 1, tmp-1, tmp-1, 1);               
              break; 
    }
    // Create Array sections
    CkMulticastMgr * mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
    uplain_proxy.ckSectionDelegate(mCastGrp);
    lplain_proxy.ckSectionDelegate(mCastGrp);
    // Create `GetBSlabsMsg' Messages 
    CkArrayIndex3D charIdx(thisIndex.x, thisIndex.y, thisIndex.z);
    CkCallback up_cb(CkIndex_JointContourNet::RecvBSlabSet1(NULL), charIdx, thisProxy);
    CkCallback lp_cb(CkIndex_JointContourNet::RecvBSlabSet2(NULL), charIdx, thisProxy);
    GetBSlabsMsg * up_msg = new GetBSlabsMsg(red_dim, true, up_cb);
    CkSetRefNum(up_msg, itr_idx);
    GetBSlabsMsg * lp_msg = new GetBSlabsMsg(red_dim, false, lp_cb);
    CkSetRefNum(lp_msg, itr_idx);
    // Multicast/Reduction on both sections 
    uplain_proxy.GetBoundSlabs(up_msg);
    lplain_proxy.GetBoundSlabs(lp_msg);

    thisProxy(thisIndex.x, thisIndex.y, thisIndex.z).Merger(itr_idx);  
  }else {
     if (lvl_id == (k >> 1)){        
       // call RecvJCN
       JCNGrMsg * gr  = new JCNGrMsg();   
       CkSetRefNum(gr, itr_idx);
       // TODO: convert local JCN to a JCNGrMsg/ Initialize gr 
       switch (red_dim){
         case 1: thisProxy( thisIndex.x - (k >> 1) , thisIndex.y, thisIndex.z).RecvJCN(gr); break; 
         case 2: thisProxy(0, thisIndex.y - (k>>1), thisIndex.z).RecvJCN(gr); break; 
         case 3: thisProxy(0, 0, thisIndex.z - (k>>1)).RecvJCN(gr); break; 
       }
//       thisProxy(thisIndex).ckDestroy();
    }
  }
}

void JointContourNet:: NextMergerStep(){
  // Set the red_dim and red_idx for next round of call(reduction)
  int k = 1 << red_idx; 
  red_idx++; 
  int cap, kk = k << 1; 
  switch (red_dim){
    case 1: cap = Xdiv; break; 
    case 2: cap = Ydiv; break;
    case 3: cap = Zdiv; break; 
  }
  if (kk > cap && !(k < cap)){
     red_idx = 1; 
     red_dim++; 
     switch (red_dim){
       case 2: if (Ydiv==1) red_dim++; 
               else { if (Zdiv == 1) fin_dim = true; 
                      break;} 
       case 3: fin_dim = true;  
               if (Zdiv==1) red_dim++; 
     }
  }
  itr_idx++; 

  if (red_dim >= 4){
    //TODO: computation ends
    finish_cb->send(); 
    // Main_Proxy.done();      
  }
}


void JointContourNet :: UpdateBorders(){
  int k = 1 << red_idx;
  // Create section on BSArray (add indices and call new)
  CProxySection_BoundSlabs sec_proxy; 
  switch (red_dim) { //TODO: check the logic one more time
    case 1: if (!fin_dim)
              sec_proxy = CProxySection_BoundSlabs::ckNew(BSArray, thisIndex.x, thisIndex.x+k-1, 1, thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1); 
            else { 
              if (thisIndex.x > 0 && thisIndex.x + k < Xdiv)
                sec_proxy = CProxySection_BoundSlabs::ckNew(BSArray, thisIndex.x, thisIndex.x+k-1, k-1, thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1); 
              else {
                if (thisIndex.x == 0)
                  sec_proxy = CProxySection_BoundSlabs::ckNew(BSArray, thisIndex.x+k-1, thisIndex.x+k-1, 1, thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1); 
                else // if (thisIndex.x + k  == Xdiv  )
                  sec_proxy = CProxySection_BoundSlabs::ckNew(BSArray, thisIndex.x, thisIndex.x, 1, thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1); 
              }
            }  
            break;
    case 2: if (!fin_dim) 
              sec_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, thisIndex.y, thisIndex.y+k-1, 1, thisIndex.z, thisIndex.z, 1); 
            else {
              if (thisIndex.y > 0 && thisIndex.y+k < Ydiv)
                sec_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, thisIndex.y, thisIndex.y+k-1, k -1, thisIndex.z, thisIndex.z, 1);
              else {
                if (thisIndex.y == 0) 
                  sec_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, thisIndex.y+k-1, thisIndex.y+k-1, 1, thisIndex.z, thisIndex.z, 1); 
                else  // if (thisIndex.y + k  == Ydiv  ) 
                  sec_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1); 
              }
            }
            break; 
    case 3: if (!fin_dim)
              sec_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, 0, Ydiv-1, 1, thisIndex.z, thisIndex.z+k-1, 1); 
            else{
              if (thisIndex.z > 0 && thisIndex.z+k < Zdiv)
                sec_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, 0, Ydiv-1, 1, thisIndex.z, thisIndex.z+k-1, k-1); 
              else{
                if (thisIndex.z == 0)    
                  sec_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, 0, Ydiv-1, 1, thisIndex.z+k-1, thisIndex.z+k-1, 1); 
                else
                  sec_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, 0, Ydiv-1, 1, thisIndex.z, thisIndex.z, 1); 
              }
            }
            break; 
  }
  CkMulticastMgr *  mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  sec_proxy.ckSectionDelegate(mCastGrp);  
  //Create  update boundary message
  CkArrayIndex3D charIdx(thisIndex.x, thisIndex.y, thisIndex.z);
  CkCallback cb(CkIndex_JointContourNet::RunMergeStep(), charIdx, thisProxy);
  UpdIDsMsg * up_msg = new (5) UpdIDsMsg(cb); // TODO: Initialize correctly
  CkSetRefNum(up_msg, itr_idx);
  // Multicast/reduction call
  sec_proxy.UpdateSlabIDs(up_msg);
}


void JointContourNet :: ComputeJCN(){}
 
void JointContourNet :: MergeJCNs(JCNGrMsg * sndGr, CkReductionMsg * fstBr, CkReductionMsg * sndBr ){
  // TODO: computation
  delete sndGr; 
  delete fstBr;
  delete sndBr;

}



void JointContourNet :: UnfoldRecvBSlabs(CkReductionMsg * msg){
  // Use CkReduction:setElement to recover each border 
  CkReduction::setElement * set_itr = (CkReduction::setElement *) msg->getData();
  while(set_itr != NULL){
    BSlabsMsg * elem = (BSlabsMsg *) &set_itr -> data; 
    // TODO : add the element to a local data structure 
    set_itr = set_itr -> next();
  }
  // TODO maybe a separate method is not necessary 
}

#include "JointContourNet.def.h"
