//
// Created by Henry Hollis on 1/28/24.
//

#ifndef LEIDEN_CCD_CCD_UTILS_H
#define LEIDEN_CCD_CCD_UTILS_H

#include <vector>
typedef std::vector<double> Vector;
class ccd_utils {
public:

    static double calcCCDsimple(const std::vector<Vector> &ref, const std::vector<Vector> &emat, bool scale);
    static std::vector<Vector> calcCorMat(const std::vector<Vector> &ref);
    static long choose(size_t n, int k);

    static void printVector(const Vector &X);

    static Vector rankVector(Vector &X);

    static float cor(const Vector &X, const Vector &Y);
};
#endif //LEIDEN_CCD_CCD_UTILS_H
