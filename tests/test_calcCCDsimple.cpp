#include "GraphHelper.h"
#include <cstdio>
#include "../include/ccd_utils.h"
using std::cout;
using std::endl;


int main() {
    //In this case, the emat is just the refCor
  double ccd = ccd_utils::calcCCDsimple(ccd_utils::refCor, 12, ccd_utils::refCor ,12,12, false );

  //Verify in deltaccd package, these are the same answers you get
  cout << "calcCCDsimple: " <<ccd <<endl;

  //Test when cor(emat) and refCor are different sizes:
    std::vector<double> emat = {
            0.154315031,	0.225517533,	0.393516789,	0.553511581,	0.352325805,
            0.433122543,	0.29048757,	0.170464952,	0.977428867	,0.311876866,
            0.59922273,	0.180537264,	0.740708598,	0.228975585,	0.186441035,
            0.995813151,	0.610829537,	0.78132706,	0.018848947	,0.920514475,
            0.826268609,	0.35945495,	0.641581163	,0.401703674	,0.430394813,
            0.650281672,	0.144778558,	0.456742456,	0.138125427,	0.158066117,
            0.242920738,	0.747656682,	0.605118059,	0.744743944,	0.381655139,
            0.456172858,	0.036161406,	0.190765146,	0.643759733,	0.840433357,
            0.060353683,	0.066428236,	0.940850344,	0.730104843,	0.123442223,
            0.82554184,	0.730472694,	0.58214786	,0.986817512,	0.87667501,
            0.094869058,	0.025067081,	0.998339339,	0.047175679,	0.366953592,
            0.239710533,	0.350125362,	0.815180148,	0.914666074,	0.86429946 //NR1D2
    };
    try{
        double ccd3 = ccd_utils::calcCCDsimple(ccd_utils::refCor,12, emat,12, 5, false);
        std::cout << "ccd3: "<<ccd3 << std::endl;
    }catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    std::vector<double> new_emat = ccd_utils::sliceColumns(emat, {0, 1, 4}, 12, 5);
    double ccd4 = ccd_utils::calcCCDsimple(ccd_utils::refCor,12, new_emat,12, 3, false);
    std::cout<<"ccd from slice cols 0, 1, 4:\n"<<ccd4<<std::endl;

    std::cout << "\nSliced Matrix:" << std::endl;
    for (size_t i = 0; i < 12; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            std::cout << new_emat[i * 3 + j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}