#ifndef COMMONDEF_H
#define COMMONDEF_H

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

#endif
