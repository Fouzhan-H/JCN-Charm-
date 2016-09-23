#ifndef COMMONDEF_H
#define COMMONDEF_H

#include "math.h"
#include <iostream>
#include <string>
#include <stdlib.h>

typedef long long int IdType; // using IdType = long; 

inline void MsgAndExit(std::string m){
   std::cout << m << std::endl; 
   exit(EXIT_FAILURE);
}  

inline void WarnMsg(std::string m){
   std::cout << m <<  std::endl; 
}

#define FP_EPSILON 1.0e-10

// Equality measure for hash map function
inline bool FP_isequal(double const& x, double const& y) 
{
  if (fabs(x - y) < FP_EPSILON)
    return 1;
  else
     return 0;
}



#endif
