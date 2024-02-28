#include "GraphHelper.h"
#include <cstdio>
#include "../include/ccd_utils.h"
using std::cout;
using std::endl;


int main() {
    // Define the dimensions of the matrix
    const size_t rows = 4;
    const size_t cols = 4;

    // Initialize the matrix with zeros
     std::vector<double> geneSampleMatrix = {1, 2, 3, 3,
                                             7, 6, 5, 7,
                                             4, 9, 8, 8,
                                             4, 0, 8, 9};

    // Print the modified matrix
    std::cout << "\nsample Matrix:" << std::endl;
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            std::cout << geneSampleMatrix[i * cols + j] << " ";
        }
        std::cout << std::endl;
    }
    std::vector<double> cormat = ccd_utils::calcCorMat(geneSampleMatrix, rows, cols);
    cout << "nrows cormat: "<< rows <<endl;
    // Print the modified matrix
    std::cout << "\nCor Matrix:" << std::endl;
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < rows; ++j) {
            std::cout << cormat[i * rows + j] << " ";
        }
        std::cout << std::endl;
    }

    std::vector<double> tmp = ccd_utils::calcCorMat(ccd_utils::refCor, 12, 12); //calculate cormat of expression matrix
    // Print the modified matrix
    std::cout << "\nrefCor cormat Matrix:" << std::endl;
    for (size_t i = 0; i < 12; ++i) {
        for (size_t j = 0; j < 12; ++j) {
            std::cout << tmp[i * 12 + j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}