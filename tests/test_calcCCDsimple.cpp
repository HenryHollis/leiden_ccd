#include "GraphHelper.h"
#include <cstdio>
#include "../include/ccd_utils.h"
using std::cout;
using std::endl;


int main() {
    const std::vector<double> refCor = {
        1.0000000,   0.77547090,  0.72492855,  0.27817942, -0.63637681, -0.60375141, -0.8614806, -0.7471112, -0.59455286, -0.8234182, -0.9146447, -0.8473980, //ARNTL
        0.7754709,   1.00000000,  0.63439613,  0.07402797, -0.62632300, -0.34987550, -0.7461844, -0.6450780, -0.70865725, -0.7845410, -0.7654845, -0.7983427, //NPAS2
        0.7249286,   0.63439613,  1.00000000,  0.06541974, -0.59727560, -0.30024636, -0.6031795, -0.6364953, -0.56958405, -0.7144612, -0.6455111, -0.7595101, //CLOCK
        0.2781794,   0.07402797,  0.06541974,  1.00000000, -0.01245765, -0.72253596, -0.4099044, -0.1411756,  0.25538496, -0.0252816, -0.3401805, -0.0781101, //CRY1
        -0.6363768, -0.62632300, -0.59727560, -0.01245765,  1.00000000,  0.28367324,  0.6234166,  0.6454257,  0.59510653,  0.6712806,  0.6618797,  0.7597038, //CRY2
        -0.6037514, -0.34987550, -0.30024636, -0.72253596,  0.28367324,  1.00000000,  0.6772739,  0.4242223, -0.06776682,  0.3366267,  0.6955807,  0.3810191, //NR1D1
        -0.8614806, -0.74618443, -0.60317949, -0.40990436,  0.62341661,  0.67727389,  1.0000000,  0.7132144,  0.52923596,  0.7673822,  0.9111478,  0.7487607, //NR1D2
        -0.7471112, -0.64507795, -0.63649530, -0.14117556,  0.64542570,  0.42422234,  0.7132144,  1.0000000,  0.60794410,  0.7467579,  0.7732704,  0.7756198, //PER1
        -0.5945529, -0.70865725, -0.56958405,  0.25538496,  0.59510653, -0.06776682,  0.5292360,  0.6079441,  1.00000000,  0.7868302,  0.5543211,  0.7530874, //PER2
        -0.8234182, -0.78454102, -0.71446119, -0.02528160,  0.67128060,  0.33662668,  0.7673822,  0.7467579,  0.78683019,  1.0000000,  0.8117621,  0.8738338, //PER3
        -0.9146447, -0.76548454, -0.64551113, -0.34018047,  0.66187971,  0.69558073,  0.9111478,  0.7732704,  0.55432112,  0.8117621,  1.0000000,  0.8443479, //DBP
        -0.8473980, -0.79834269, -0.75951011, -0.07811010,  0.75970381,  0.38101906,  0.7487607,  0.7756198,  0.75308740,  0.8738338,  0.8443479,  1.0000000 //TEF
};
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
        double ccd = ccd_utils::calcCCDsimple(refCor,12, emat,12, 5, false);
        std::cout << "ccd: "<<ccd << std::endl;
    }catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    std::vector<double> new_emat = ccd_utils::sliceColumns(emat, {0, 1, 4,6}, 12, 5);
    double ccd2 = ccd_utils::calcCCDsimple(refCor,12, new_emat,12, 3, false);
    std::cout<<"ccd from slice cols 0, 1, 4:\n"<<ccd2<<std::endl;

    std::vector<double> new_emat2 = ccd_utils::sliceColumns(emat, {0, 4, 1}, 12, 5);
    double ccd3 = ccd_utils::calcCCDsimple(refCor,12, new_emat2,12, 3, false);
    std::cout<<"ccd from slice cols 0, 4, 1:\n"<<ccd3<<std::endl;

//    std::cout << "\nSliced Matrix:" << std::endl;
//    for (size_t i = 0; i < 12; ++i) {
//        for (size_t j = 0; j < 3; ++j) {
//            std::cout << new_emat[i * 3 + j] << " ";
//        }
//        std::cout << std::endl;
//    }

    return 0;
}