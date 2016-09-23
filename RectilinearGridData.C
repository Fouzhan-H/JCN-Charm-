#include "RectilinearGridData.h"
#include <stdio.h>
#include <iostream>
#include <fstream>

RectilinearGridData::RectilinearGridData(int M, int * start, int* size, const int* glb_dim, int N, char fileNames [][100])
  : st(std::vector<IdType>(M))
  , sz(std::vector<int>(M)), 
  M_geom (M), N(N), glb_dim(glb_dim), fnames(fileNames){

  for (int i = 0; i < M_geom; i++){
    st[i] = *(start+i);
    sz[i] = *(size+i);
  }

  cCell = -1;       
}


void RectilinearGridData::fetchInput(int i, float * mem_block){
  const int file_flt_sz = 4;   // Each sample is stored as a 4-byte floating point
  std::streampos gpos, gposl; 
  // The dataset dimension 
  IdType gx_dim = * glb_dim; 
  IdType gy_dim = * (glb_dim + 1); 
  IdType lx_st, ly_st, lz_st; 
  int lx_dim, ly_dim, lz_dim; 

  //set the local vars for 2 and 3 dimensions  
  lx_st = st[0]; lx_dim = sz[0]; 
  ly_st = st[1]; ly_dim = sz[1]; 
  if (M_geom == 2){
     lz_st = 0; 
     lz_dim = 1;
  }
  else {
    lz_st = st[2]; 
    lz_dim = sz[2]; 
  }


  // Input files store the whold dataset,
  // To read a subset of data, first, the index of 
  // first element of the subset mush be computed.   
  gpos = (lz_st * gx_dim * gy_dim + 
          ly_st * gx_dim + 
          lx_st ) * file_flt_sz;
   
  std::ifstream inf (fnames[i], std::ios::binary);    // open input file
    
  // readdata
  float * mem_block_ptr = mem_block; 
  // Samples in the files are arranged in increasing 
  // x-dimension, then along y and z dimensions,  
  // respectively. So data is read using two nested 
  // for loop on Z and Y dimension. The read index 
  // should be moved to the correct position after 
  // each iteration.
  
  for (int i = 0; i < lz_dim; i++){
    inf.seekg(gpos);         // move the file read index to the start of the plane 
    gposl = gpos; 
    for (int j = 0; j < ly_dim; j++){
      inf.seekg(gposl);
      inf.read(reinterpret_cast<char *> (mem_block_ptr), lx_dim*file_flt_sz); // read one line in x dimension 
      mem_block_ptr += lx_dim ;
      gposl += file_flt_sz * gx_dim; 
    } //y-dimension
    gpos += file_flt_sz * gx_dim * gy_dim; 
  } //z-dimension

  inf.close(); // close the input file
  
}


IdType RectilinearGridData::GetCellsNr(){
  IdType cnr = 1; 
  for (int i = 0; i < M_geom; i++)
       cnr *= sz[i] - 1; 
  
  if (M_geom == 2) 
    return 2 * cnr; 

  return 6*cnr; 
} 


int RectilinearGridData::CellComponentNr(){
  if (M_geom == 2) 
     return 3; 
  return 4; 
}

void RectilinearGridData::GetSpatialBounds(double * bnds){
  for (int i = 0; i < M_geom; i++){
     bnds[2*i] = st[i]; 
     bnds[2*i+1] = st[i] + sz[i] - 1;  
  }
}

IdType RectilinearGridData::GetCellComponent(IdType c, int j , double * point){
  int cx, cy, cz;
  int cr, cyr, cs; 
  IdType ridx; 

  // Compute the index of cubic cell which 
  // include cell "c". These values will be 
  // cached az cell_x, cell_y, and cell_z. 
  if (M_geom == 2){
    if (c != cCell){
      int xmax = sz[0] - 1;
      int cellPerRow = 2 * xmax; 

      cell_z = 0;
      cell_y = c / cellPerRow; 
      cyr = c % cellPerRow; 
      cell_x = cyr / 2; 
      cs = cyr % 2;
      Simplicate3D(cs);    
    }
  } else {
     if (c != cCell){
       int xmax = sz[0] - 1;
       int ymax = sz[1] - 1; 
       int cellPerPlane = 6 * xmax * ymax;     
       int cellPerRow = 6 * xmax; 

       cell_z = c / cellPerPlane;
       cr = c % cellPerPlane;
       cell_y = cr / cellPerRow; 
       cyr = cr % cellPerRow; 
       cell_x = cyr / 6; 
       cs = cyr % 6;
       Simplicate3D(cs);    
     }
  }
  
  cx = cell_x;
  cy = cell_y;
  cz = cell_z;
  cCell = c; 

  switch (cell4[j]){
    case 0: break ; 
    case 1: cx++; 
            break;
    case 2: cy++; 
            break;   
    case 3: cx++; cy++; 
            break;
    case 4: cz++; 
            break;  
    case 5: cz++; cx++ ; 
            break;
    case 6: cz++; cy++; 
            break;   
    case 7: cz++; cx++; cy++; 
            break;
  }

  point[0] = st[0] + cx;
  point[1] = st[1] + cy;
  point[2] = st[2] + cz;
 
 
  ridx = cx + cy*sz[0] + cz*sz[0]*sz[1];
  return ridx; 
}


void RectilinearGridData::Simplicate3D(short cs ){
  if (M_geom == 2){ 
     switch (cs){
       case 0: cell4[0] = 0; cell4[1] = 1; cell4[2]= 3; cell4[3]=0;  // cell4 = {0, 1, 3, 0};
               break; 
       case 1: cell4[0] = 0; cell4[1] = 2; cell4[2]= 3; cell4[3]=0;  // cell4 = {0, 2, 3, 0};  
               break;
     }
     return; 
  }
   
  switch (cs){
    case 0:  cell4[0] = 2; cell4[1] = 3; cell4[2]= 0; cell4[3]=7;    // cell4 = {2, 3, 0, 7};
            break;
    case 1:  cell4[0] = 2; cell4[1] = 0; cell4[2]= 6; cell4[3]=7;    // cell4 = {2, 0, 6, 7};
            break; 
    case 2:  cell4[0] = 3; cell4[1] = 0; cell4[2]= 1; cell4[3]=7;    // cell4 = {3, 0, 1, 7};
            break; 
    case 3:  cell4[0] = 0; cell4[1] = 1; cell4[2]= 7; cell4[3]=5;    // cell4 = {0, 1, 7, 5};
            break; 
    case 4:  cell4[0] = 0; cell4[1] = 6; cell4[2]= 7; cell4[3]=4;    // cell4 = {0, 6, 7, 4};
            break; 
    case 5:  cell4[0] = 0; cell4[1] = 7; cell4[2]= 4; cell4[3]=5;    // cell4 = {0, 7, 4, 5};
            break; 
  }
  
}



void RectilinearGridData::fetchCharInput(int i, float * mem_block){
  const int file_flt_sz = 1;   // Each sample is stored as a 4-byte floating point
  std::streampos gpos, gposl; 
  // The dataset dimension 
  IdType gx_dim = * glb_dim; 
  IdType gy_dim = * (glb_dim + 1); 
  IdType lx_st, ly_st, lz_st; 
  int lx_dim, ly_dim, lz_dim; 

  //set the local vars for 2 and 3 dimensions  
  lx_st = st[0]; lx_dim = sz[0]; 
  ly_st = st[1]; ly_dim = sz[1]; 
  if (M_geom == 2){
     lz_st = 0; 
     lz_dim = 1;
  }
  else {
    lz_st = st[2]; 
    lz_dim = sz[2]; 
  }


  // Input files store the whold dataset,
  // To read a subset of data, first, the index of 
  // first element of the subset mush be computed.   
  gpos = (lz_st * gx_dim * gy_dim + 
          ly_st * gx_dim + 
          lx_st ) * file_flt_sz;
   
  std::ifstream inf (fnames[i], std::ios::binary);    // open input file
    
  // readdata
  char * mem_block_ptr = new char [lx_dim]; 
  float *fmem_block_ptr = mem_block; 
  // Samples in the files are arranged in increasing 
  // x-dimension, then along y and z dimensions,  
  // respectively. So data is read using two nested 
  // for loop on Z and Y dimension. The read index 
  // should be moved to the correct position after 
  // each iteration.
  
  for (int i = 0; i < lz_dim; i++){
    inf.seekg(gpos);         // move the file read index to the start of the plane 
    gposl = gpos; 
    for (int j = 0; j < ly_dim; j++){
      inf.seekg(gposl);
      inf.read(mem_block_ptr, lx_dim*file_flt_sz); // read one line in x dimension 
      for (int k = 0; k < lx_dim; k++ )
          fmem_block_ptr[k] = (unsigned char) mem_block_ptr[k]; 
      fmem_block_ptr += lx_dim ;
      gposl += file_flt_sz * gx_dim; 
    } //y-dimension
    gpos += file_flt_sz * gx_dim * gy_dim; 
  } //z-dimension

  inf.close(); // close the input file
  delete mem_block_ptr; 
}

