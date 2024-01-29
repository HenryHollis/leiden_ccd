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

  //Test when cor(emat) and refCor are different sizes:
    std::vector<Vector> emat = {
            {1.0000000,   0.77547090,  0.72492855,  0.27817942, -0.63637681, -0.60375141, -0.8614806, -0.7471112, -0.59455286, -0.8234182, -0.9146447, -0.8473980}, //ARNTL
            {0.7754709,   1.00000000,  0.63439613,  0.07402797, -0.62632300, -0.34987550, -0.7461844, -0.6450780, -0.70865725, -0.7845410, -0.7654845, -0.7983427}, //NPAS2
            {0.7249286,   0.63439613,  1.00000000,  0.06541974, -0.59727560, -0.30024636, -0.6031795, -0.6364953, -0.56958405, -0.7144612, -0.6455111, -0.7595101}, //CLOCK
            {0.2781794,   0.07402797,  0.06541974,  1.00000000, -0.01245765, -0.72253596, -0.4099044, -0.1411756,  0.25538496, -0.0252816, -0.3401805, -0.0781101}, //CRY1
            {-0.6363768, -0.62632300, -0.59727560, -0.01245765,  1.00000000,  0.28367324,  0.6234166,  0.6454257,  0.59510653,  0.6712806,  0.6618797,  0.7597038}, //CRY2
            {-0.6037514, -0.34987550, -0.30024636, -0.72253596,  0.28367324,  1.00000000,  0.6772739,  0.4242223, -0.06776682,  0.3366267,  0.6955807,  0.3810191}, //NR1D1
            {-0.8614806, -0.74618443, -0.60317949, -0.40990436,  0.62341661,  0.67727389,  1.0000000,  0.7132144,  0.52923596,  0.7673822,  0.9111478,  0.7487607}, //NR1D2
    };
    try{
        double ccd3 = ccd_utils::calcCCDsimple(ccd_utils::refCor,emat, false);

    }catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}