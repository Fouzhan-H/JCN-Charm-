#include "Main.decl.h"
#include "Main.h"

#include "JointContourNet.h"
#include "BoundSlabs.h"

#include <ckmulticast.h>

Main :: Main(CkArgMsg * msg){
  startTime = CkWallTimer(); 
 
  GetParams(msg); // Read Program Arguments: 
  
  mCastGrpId = CProxy_CkMulticastMgr::ckNew();
  // Create JCN and BoundSlab Arrays
  CProxy_JointContourNet JCNArray = CProxy_JointContourNet::ckNew( Xdiv, Ydiv, Zdiv);
  CkArrayOptions opts(Xdiv, Ydiv, Zdiv); 
  opts.bindTo(JCNArray); 
  BSArray = CProxy_BoundSlabs::ckNew(3, opts);
  JCNArray.Start(); 



/*  CkArrayOptions opts(Xdiv, Ydiv, Zdiv); 
  BSArray = CProxy_BoundSlabs::ckNew(3, opts);
  CProxy_JointContourNet JCNArray = CProxy_JointContourNet::ckNew( Xdiv, Ydiv, Zdiv);
  JCNArray.Start(); */
}

/*void Main :: Done(){
  double endTime = CkWallTimer();
  CkPrintf("all done in %f seconds.\n", endTime - startTime);
  CkExit();
}*/

// Read Program Arguments: 
//   * Dataset dimension 
//   * Slab widths 
//   * Slab bases
//   * division across dimensions 
void Main :: GetParams (CkArgMsg * m){
  bool dimp = false, divp = false, rngp = false, swp = false, bsp = false, fnp=false; 

  for (int i = 1 ; i < m -> argc; ){
     if (strcmp("-d", m->argv[i]) == 0){ //Read and set dataset dimension
        XDim = atoi(m->argv[++i]);
        YDim = atoi(m->argv[++i]);
        ZDim = atoi(m->argv[++i]);
        i++; dimp = true; 
        continue;
     }

     if (strcmp("-rd", m->argv[i]) == 0){ // Read and set the range domain dimension
       Range_Dim = atoi(m -> argv[++i]);
       i++; rngp = true;   
       continue; 
     }

     if (strcmp("-sw", m->argv[i]) == 0){ //Read and set slab width for each range dimension
       if (!rngp) {
          CkPrintf("Range dimension (-rd) is expected before slab width values\n");
          CkExit();
       }
       for (int j = 0; j < Range_Dim; j++){  
         sws[j] = atof(m->argv[++i]); 
       } 
       i++; swp = true; 
       continue; 
     }

     if (strcmp("-bs", m->argv[i]) == 0){ // Read and set Slab bases for each dimension 
       if (!rngp) {
          CkPrintf("Range dimension (-rd) is expected before slab bases\n");
          CkExit();
       }
       for (int j = 0 ; j < Range_Dim; j++)
         bss[j] = atof(m->argv[++i]);
       i++; bsp = true; 
       continue; 
     } 
     
     if (strcmp("-f", m->argv[i]) == 0 ){
       if (!rngp) {
          CkPrintf("Range dimension (-rd) is expected before list of input file names\n");
          CkExit();
       }
       for (int j = 0 ; j < Range_Dim; j++)
         strcpy(fns[j],m->argv[++i]);
       i++; fnp = true; 
       continue; 
     
     }
 
     if (strcmp("-div", m->argv[i]) == 0){ // Read and set number of division for each domain dimension  
       Xdiv = atoi(m->argv[++i]);
       Ydiv = atoi(m->argv[++i]);
       Zdiv = atoi(m->argv[++i]);
       i++;  divp = true; 
       continue; 
     } 

     break;
  }

  if (!dimp) CkPrintf("Dataset dimension (-d) is expected \n");
  if (!rngp) CkPrintf("Range dimension (-rd) is expected \n");
  if (!swp) CkPrintf("Slab widths (-sw) are expected \n");
  if (!fnp) CkPrintf("Input file names (dataset) is expected\n"); 
 
  if (!dimp || !rngp || !swp || !fnp) CkExit(); 

  if (!bsp){ // Set default value for slab bases to zero
    CkPrintf ("Warning: Slabs are computed with base 0 \n"); 
    for (int j =0; j < Range_Dim; j++)
      bss[j] = 0; 
  }

  if (!divp){ // Set default value for domain division to one, no parallelism
    CkPrintf ("Warning: Compution is done with no parallelism \n"); 
    Xdiv = Ydiv = Zdiv = 1; 
  }
  
}

Main :: Main(CkMigrateMessage * msg){}

#include "Main.def.h"

