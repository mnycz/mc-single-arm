≈#ifdef __cplusplus
extern "C"
#endif

#include <utility> //std::pair                                                                                                                                                                          
#include <iostream>
#include <fstream>
#include <iomanip>
//#include <TROOT.h>
//#include <TMath.h>




double Landau(double x, double y){
  double Random_Number=0.0;
  //gRandom = new TRandom3;
 
  Random_Number =gRandom-> Landau(x,y);
  return Random_Number;
}
