#include "JointContourNet.decl.h"
#include "JointContourNet.h"

#include "RectilinearGridData.h"
#include "CompJointContourNet.h"
#include  "MergJointContourNet.h"

#include <fstream>

JointContourNet :: JointContourNet(): finish_cb(NULL), slabs(NULL), edges(NULL){
  fin_dim = (Ydiv == 1 && Zdiv == 1) ? true: false; 
  itr_idx = 0;  
  red_dim =1 ; red_idx = 0; 
  NextMergerStep(); 
}
 
JointContourNet :: JointContourNet(CkMigrateMessage * msg){
  finish_cb = NULL; 
  slabs = NULL;
  edges = NULL; 
}

JointContourNet :: ~JointContourNet(){
  if (finish_cb != NULL)
    delete finish_cb;
  if (slabs != NULL)
    delete slabs; 
  if (edges != NULL)
    delete edges; 
}

void JointContourNet::pup(PUP::er &p){
  CBase_JointContourNet::pup(p);
  __sdag_pup(p);

  p| itr_idx;
  p| red_dim;
  p| red_idx; 
  p| fin_dim; 
  p| upd_sec1;
  p| upd_sec2;
  //PUP data member finish_cb of type CkCallback* 
  int has_f_cb = (finish_cb != NULL);
  p| has_f_cb; 
  if (has_f_cb){
    if (p.isUnpacking()) finish_cb = new CkCallback();
    p | *finish_cb; 
  }else 
     finish_cb = NULL; 
  //PUP data member slabs; 
  int has_slabs = (slabs != NULL);
  p| has_slabs; 
  if (has_slabs){
    if (p.isUnpacking()) slabs = new std::vector<std::vector<double>> (0, std::vector<double>(Range_Dim));
    p | *slabs; 
  }else 
     slabs = NULL;
  //PUUP data member edges
  int has_edges = (edges != NULL);
  p| has_edges; 
  if (has_edges){
    if (p.isUnpacking()) edges = new std::set<std::pair<long long, long long>>();
    p| *edges; 
  }else 
     edges = NULL; 
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

  int gz[3] = {XDim, YDim, ZDim};

  int M_geom = 3; //TODO 
  RectilinearGridData dset (M_geom, st, sz, gz , Range_Dim , fns); //create RecGrid  
  //TODO: call ComputeJCN() 
  long long slabNr, bndSlabNr;  
  long long pointNr = dim_x * dim_y * dim_z; 
  CompJointContourNet cmpJCN (&dset, Range_Dim, bss, sws, pointNr); //TODO

  cmpJCN.compJCN(bndSlabNr, slabNr); 
 
  slabs = new std::vector<std::vector<double>> (slabNr, std::vector<double>(Range_Dim));
  edges = new std::set<std::pair<long long, long long>>();
  //create populate message 
  CkCallback mrg_cb(CkIndex_JointContourNet::RunMergeStep()
                   , CkArrayIndex3D(thisIndex.x, thisIndex.y, thisIndex.z), thisProxy);
  BSPopulateMsg * msg = new (M_geom, M_geom, bndSlabNr*M_geom, bndSlabNr ) BSPopulateMsg(mrg_cb, slabNr, bndSlabNr);
  for (int i = 0; i < M_geom; i++){
     msg -> st [i] = st[i]; 
     msg -> sz [i] = sz[i]; 
  }
  cmpJCN.extractJCN(*slabs, *edges, msg->facetCntrs, msg->slabIDs); 


  //TODO Populate the coressponding BoundSlab Array element   
  BoundSlabs * bslocal = BSArray(thisIndex.x, thisIndex.y, thisIndex.z).ckLocal();
 if (bslocal == NULL)
     BSArray(thisIndex.x, thisIndex.y, thisIndex.z).Populate(msg);
  else  
    bslocal->Populate(msg);

}  

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
    // Create `SndBndMsg' Messages 
    CkArrayIndex3D charIdx(thisIndex.x, thisIndex.y, thisIndex.z);
    CkCallback up_cb(CkIndex_JointContourNet::RecvBSlabSet1(NULL), charIdx, thisProxy);
//    CkCallback lp_cb(CkIndex_JointContourNet::RecvBSlabSet2(NULL), charIdx, thisProxy);
    SndBndMsg * up_msg = new SndBndMsg(red_dim, true, up_cb);
    CkSetRefNum(up_msg, itr_idx);
    TargetBndMsg * lp_msg = new TargetBndMsg(red_dim, false);
    CkSetRefNum(lp_msg, itr_idx);
    // Multicast/Reduction on both sections 
    uplain_proxy.GetAdjacentSlabs(up_msg);
    lplain_proxy.SndBndFaces(lp_msg);
    thisProxy(thisIndex.x, thisIndex.y, thisIndex.z).Merger(itr_idx);  

  }else {
     if (lvl_id == (k >> 1)){        
       // Convert local JCN to a JCNGrMsg/ Initialize gr 
       int slabsSz = Range_Dim * slabs->size();
       int edgesSz = 2 * edges->size(); 
       JCNGrMsg * gr  = new (slabsSz, edgesSz ) JCNGrMsg();   
       // Copy slab values into the JCNGrMsg
       double * gr_vexp = gr -> vexs; 
       for (int i = 0; i < slabs -> size(); i++){
         for (short j = 0 ; j < Range_Dim; j++)
           *gr_vexp++= (*slabs)[i][j];
       }    
       // Copy edges into the JCNGrMsg
       int edge_it = 0; 
       for (std::set<std::pair<long long, long long>>::iterator it = edges->begin(); it != edges->end(); ++it, edge_it+=2 ){
          gr->edgs [edge_it] = std::get<0>(*it); 
          gr->edgs [edge_it +1] = std::get<1>(*it); 
       }
       gr -> slabCnt = slabs->size();
       gr -> edgCnt = edges->size();
       CkSetRefNum(gr, itr_idx);           //set the message tag (reduction Level)
       switch (red_dim){
         case 1: thisProxy( thisIndex.x - (k >> 1) , thisIndex.y, thisIndex.z).RecvJCN(gr); break; 
         case 2: thisProxy(0, thisIndex.y - (k>>1), thisIndex.z).RecvJCN(gr); break; 
         case 3: thisProxy(0, 0, thisIndex.z - (k>>1)).RecvJCN(gr); break; 
       }
       delete slabs;
       slabs = NULL; 
       delete edges; 
       edges = NULL; 
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
    // write out final JCN in a file
    CkPrintf("Slab Nr: %i Edge Nr: %i \n", slabs -> size(), edges -> size());
    JCNtoFile("jcn.vtk");
    finish_cb->send(); 
    // Main_Proxy.done();      
  }
}


void JointContourNet :: UpdateBorders(int set1_Nr, long long * set1_ids, int set2_Nr, long long * set2_ids){
  int k = 1 << red_idx;
  int hk = 1 << (red_idx - 1);
  upd_sec1 = true; 
  upd_sec2 = true; 
  // Create section on BSArray (add indices and call new)
  CProxySection_BoundSlabs sec1_proxy; 
  CProxySection_BoundSlabs sec2_proxy; 
  switch (red_dim) { //TODO: check the logic one more time
    case 1: if (!fin_dim){
              sec1_proxy = CProxySection_BoundSlabs::ckNew(BSArray, thisIndex.x, thisIndex.x + hk-1, 1, thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1); 
              sec2_proxy = CProxySection_BoundSlabs::ckNew(BSArray, thisIndex.x + hk, thisIndex.x+k-1, 1, thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1); 
            }else { 
               if (thisIndex.x > 0 && thisIndex.x + k < Xdiv){
                 sec1_proxy = CProxySection_BoundSlabs::ckNew(BSArray, thisIndex.x, thisIndex.x, 1, thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1); 
                 sec2_proxy = CProxySection_BoundSlabs::ckNew(BSArray, thisIndex.x+k-1, thisIndex.x+k-1, 1, thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1); 
               } else {
                  if (thisIndex.x == 0){
                    upd_sec1= false; 
                    sec2_proxy = CProxySection_BoundSlabs::ckNew(BSArray, thisIndex.x+k-1, thisIndex.x+k-1, 1, thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1); 
                  }else{ // if (thisIndex.x + k  == Xdiv  )
                    sec1_proxy = CProxySection_BoundSlabs::ckNew(BSArray, thisIndex.x, thisIndex.x, 1, thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1); 
                    upd_sec2 = false; 
                  }
              }
            }  
            break;
    case 2: if (!fin_dim){ 
              sec1_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, thisIndex.y, thisIndex.y+hk-1, 1, thisIndex.z, thisIndex.z, 1); 
              sec2_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, thisIndex.y+hk, thisIndex.y+k-1, 1, thisIndex.z, thisIndex.z, 1); 
            } else {
               if (thisIndex.y > 0 && thisIndex.y+k < Ydiv){
                sec1_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1);
                sec2_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, thisIndex.y+k-1, thisIndex.y+k-1, 1, thisIndex.z, thisIndex.z, 1);
               } else {
                  if (thisIndex.y == 0){ 
                    upd_sec1 = false; 
                    sec2_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, thisIndex.y+k-1, thisIndex.y+k-1, 1, thisIndex.z, thisIndex.z, 1); 
                  }else{  // if (thisIndex.y + k  == Ydiv  ) 
                    sec1_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, thisIndex.y, thisIndex.y, 1, thisIndex.z, thisIndex.z, 1); 
                    upd_sec2 = false; 
                  }
               }
            }
            break; 
    case 3: if (!fin_dim){
              sec1_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, 0, Ydiv-1, 1, thisIndex.z, thisIndex.z+hk-1, 1); 
              sec2_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, 0, Ydiv-1, 1, thisIndex.z+hk, thisIndex.z+k-1, 1); 
            } else{
              if (thisIndex.z > 0 && thisIndex.z+k < Zdiv){
                sec1_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, 0, Ydiv-1, 1, thisIndex.z, thisIndex.z, 1); 
                sec2_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, 0, Ydiv-1, 1, thisIndex.z+k-1, thisIndex.z+k-1, 1); 
              } else{
                if (thisIndex.z == 0){    
                  upd_sec1 = false; 
                  sec2_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, 0, Ydiv-1, 1, thisIndex.z+k-1, thisIndex.z+k-1, 1); 
                }else{
                  sec1_proxy = CProxySection_BoundSlabs::ckNew(BSArray, 0, Xdiv-1, 1, 0, Ydiv-1, 1, thisIndex.z, thisIndex.z, 1); 
                  upd_sec2 = false; 
                }
              }
            }
            break; 
  }
  CkMulticastMgr *  mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  //Create  callback
  CkArrayIndex3D charIdx(thisIndex.x, thisIndex.y, thisIndex.z);

  long long * srcp, *desp; 
  if (upd_sec1){
    sec1_proxy.ckSectionDelegate(mCastGrp);  
    CkCallback cb(CkIndex_JointContourNet::Updated1(), charIdx, thisProxy);
    UpdIDsMsg * up_msg1 = new (set1_Nr) UpdIDsMsg(cb);
    std::copy(set1_ids, set1_ids+set1_Nr, up_msg1 -> ids); // Copy id mapping from Union-Find structure into an UpdIDsMsg 
    CkSetRefNum(up_msg1, itr_idx);
    sec1_proxy.UpdateSlabIDs(up_msg1);                     // Multicast/reduction call
  }
  
  if (upd_sec2){
    sec2_proxy.ckSectionDelegate(mCastGrp);  
    CkCallback cb(CkIndex_JointContourNet::Updated2(), charIdx, thisProxy);
    UpdIDsMsg * up_msg2 = new (set2_Nr) UpdIDsMsg(cb);
    std::copy(set2_ids, set2_ids + set2_Nr, up_msg2->ids); // Copy id mapping from Union-Find structure into an UpdIDsMsg 
    CkSetRefNum(up_msg2, itr_idx);
    sec2_proxy.UpdateSlabIDs(up_msg2);                     // Multicast/reduction call
  }
}

void JointContourNet :: ComputeJCN(){}
 
void JointContourNet :: MergeJCNs(JCNGrMsg * rcvdGr, CkReductionMsg * adjSlabs){
  
  std::vector<std::vector<double>> * slabs_gr1 = slabs;
  std::set<std::pair<long long, long long>> * edges_gr1 = edges;

  MergJointContourNet jcn_merger (slabs_gr1 -> size(), slabs_gr1, rcvdGr -> slabCnt, rcvdGr -> vexs, Range_Dim); 
  
  // Use CkReduction:setElement to recover each shared border 
  // And process adjacent slabs by calling AddAdjacentSlabs(..) of jcn_merger
  CkReduction::setElement * set_itr = (CkReduction::setElement *) adjSlabs->getData();
  long long * elem;
  while(set_itr != NULL){
    elem = (long long *) &set_itr -> data; 
    jcn_merger.AddAdjacentSlabs(elem + 1 , *elem);     
    set_itr = set_itr -> next();
  }

  long long slabsNr;
  jcn_merger.ExtractJCNSlabs(slabsNr);

  slabs = new std::vector<std::vector<double>>(slabsNr, std::vector<double>(Range_Dim)); 
  edges = new std::set<std::pair<long long, long long>>();

  jcn_merger.ExtractJCNGr(*slabs, *edges, *edges_gr1, rcvdGr -> edgs, rcvdGr -> edgCnt);  

  long long * uids = jcn_merger.UpdIdMap();
  UpdateBorders(slabs_gr1 -> size(), uids, rcvdGr -> slabCnt, uids + (slabs_gr1 -> size()));

  delete slabs_gr1; 
  delete edges_gr1; 
  delete rcvdGr; 
  delete adjSlabs;
}

void JointContourNet::JCNtoFile(std::string fname){
   std::ofstream fh (fname);

   IdType vn = slabs->size();
   IdType en = edges->size(); 
   int rngn = slabs->operator[](0).size();

   fh << "# vtk DataFile Version 3.0" << std::endl 
      << "vtk output" << std::endl
      << "ASCII" << std::endl
      << "DATASET UNDIRECTED_GRAPH"  << std::endl
      << "POINTS " << vn << " float" << std::endl; 
      

   // print out the domain coordinates 
   IdType td = vn / 3; 
   int tr = vn % 3; 
   for (IdType i = 0; i < td ; i++)
       fh << "0 0 0 0 0 0 0 0 0" << std::endl; 
   if (tr != 0) {
     for (int i = 0; i < tr; i++)
         fh << "0 0 0 "; 
     fh << std::endl; 
   }
   
   fh << "VERTICES " <<  vn << std::endl; 
   fh << "EDGES " << en << std::endl; 

   // print out edges. 
   for (std::set<std::pair<long long, long long>>::iterator it=edges->begin(); it != edges ->end(); ++it)
     fh << (*it).first << " " << (*it).second << std::endl; 

   fh << "VERTEX_DATA " << vn << std::endl; 
   fh << "FIELD FieldData " <<  rngn << std::endl; 


   // print out the slab values    
   for (int i =0; i < rngn; i++){
     fh << fns[i] << " 1 " << vn << " double"   << std::endl;
 
     for (IdType j = 0; j < vn;){
        for (IdType  k = 0; j < vn && k < 10 ; j++, k++ ){           
           fh << (*slabs)[j][i] << " "; 
        }
        fh << std::endl; 
     }
   } // loop on range dimension
   
  fh.close();
}

#include "JointContourNet.def.h"
