#include "GraphHelper.h"
#include <cstdio>
#include "../include/ccd_utils.h"
using std::cout;
using std::endl;


int main() {
    //In this case, the emat is just the refCor
  double ccd = ccd_utils::calcCCDsimple(ccd_utils::refCor, ccd_utils::refCor , false);

//Also test the slice columns functionality, slicing the first 6 columns of the refcor to form a new emat
  double ccd2 = ccd_utils::calcCCDsimple(ccd_utils::refCor, ccd_utils::sliceColumns(ccd_utils::refCor, {0,1,2,3,4,5}), false);

  //Verify in deltaccd package, these are the same answers you get
  cout << "calcCCDsimple: " <<ccd <<endl;
  cout << "ccd of first 6 cols of emat: " <<ccd2 <<endl;

  return 0;
}