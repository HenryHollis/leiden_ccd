#include "GraphHelper.h"
#include <cstdio>
#include "../include/ccd_utils.h"
using std::cout;
using std::endl;


int main() {
    // Define the dimensions of the matrix
    const size_t rows = 4;
    const size_t cols = 3;

    // Initialize the matrix with zeros
    std::vector<Vector> geneSampleMatrix(rows, Vector(cols, 0.0));

    geneSampleMatrix[0] = {-1.22823867, -0.096354563,-0.4288392};
    geneSampleMatrix[1] = {0.04329142,1.390203721,0.0412204};
    geneSampleMatrix[2] = {0.69389058, 0.007663258, 0.4592899};
    geneSampleMatrix[3] = {-0.68343883,-0.410407624, -0.2839379};

    // Print the modified matrix
    std::cout << "\nsample Matrix:" << std::endl;
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            std::cout << geneSampleMatrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    vector<Vector> cormat = ccd_utils::calcCorMat(geneSampleMatrix);
    cout << "nrows cormat: "<< cormat.size() <<endl;
    // Print the modified matrix
    std::cout << "\nCor Matrix:" << std::endl;
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < rows; ++j) {
            std::cout << cormat[i][j] << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}