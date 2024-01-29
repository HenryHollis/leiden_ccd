//
// Created by Henry Hollis on 1/28/24.
//

#include "../include/ccd_utils.h"
#include <iostream>
#include <cmath>
#include <stdexcept>


double ccd_utils::calcCCDsimple(const std::vector<Vector> &ref,
                     const std::vector<Vector> &emat,
                      bool scale) {

    std::vector<Vector> cormat = calcCorMat(emat); //calculate cormat of expression matrix
    if (cormat.size() != ref.size() || cormat.empty() || ref[0].size() != cormat[0].size()) {
        throw std::invalid_argument("Matrices must be of the same size for calcCCDsimple");
    }
    size_t numRows = ref.size();
    double upperTriDiff = 0.0;
    //loop through triangular matrix and accumulate the difference between entries of ref and cormat
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = i; j < numRows; ++j) {
            upperTriDiff += pow(ref[i][j] - cormat[i][j], 2);
        }
    }
    double ccd = sqrt(upperTriDiff);

    if (!ref.empty() && !ref[0].empty()) {
        // Number of columns is the size of any row (assuming all rows have the same size)
        size_t numColumns = ref[0].size();
        if (scale) {
            size_t nPairs = choose(numColumns, 2);
            ccd /= static_cast<double>(nPairs);
        }
    } else
        std::cerr << "Matrix ref is empty or has empty rows." << std::endl;

    return ccd;
}

std::vector<Vector> ccd_utils::calcCorMat(const std::vector<Vector> &ref) {
   // size_t numCols = ref[0].size();
    size_t numRows = ref.size();

    // Initialize the correlation matrix with zeros
    std::vector<Vector> correlationMatrix(numRows, Vector(numRows, 0.0));

    // Calculate pairwise correlations
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = i; j < numRows; ++j) {
            if (i == j) {
                // Diagonal elements (correlation with itself) are always 1
                correlationMatrix[i][j] = 1.0;
            } else {
                // Off-diagonal elements are calculated using the cor function
                Vector rank_x = rankVector(const_cast<Vector &>(ref[i])); //rank vectors first to calculate spearman coeff
                Vector rank_y = rankVector(const_cast<Vector &>(ref[j]));

                correlationMatrix[i][j] = correlationMatrix[j][i] = cor(rank_x, rank_y);
            }
        }
    }

    return correlationMatrix;
}


long ccd_utils::choose(size_t n, int k) {
    if (0 == k)
        return 1;
    return (n * choose(n - 1, k - 1)) / k;
}

/**
 * Spearman Correlation Code
 * Inspired from https://www.geeksforgeeks.org/program-spearmans-rank-correlation/
 */
// Utility Function to print vector
void ccd_utils::printVector(const Vector &X) {
    for (auto i: X)
        std::cout <<i<<" ";
    std::cout << std::endl;
}

// Function returns the rank vector
// of the set of observations
Vector ccd_utils::rankVector(Vector &X) {

    int N = X.size();

    // Rank Vector
    Vector Rank_X(N);

    for(int i = 0; i < N; i++)
    {
        int r = 1, s = 1;

        // Count no of smaller elements
        // in 0 to i-1
        for(int j = 0; j < i; j++) {
            if (X[j] < X[i] ) r++;
            if (X[j] == X[i] ) s++;
        }

        // Count no of smaller elements
        // in i+1 to N-1
        for (int j = i+1; j < N; j++) {
            if (X[j] < X[i] ) r++;
            if (X[j] == X[i] ) s++;
        }

        // Use Fractional Rank formula
        // fractional_rank = r + (n-1)/2
        Rank_X[i] = r + (s-1) * 0.5;
    }

    // Return Rank Vector
    return Rank_X;
}

// function that returns
// Pearson correlation coefficient.
float ccd_utils::cor(const Vector &X, const Vector &Y) {
    int n = X.size();
    float sum_X = 0, sum_Y = 0,
            sum_XY = 0;
    float squareSum_X = 0,
            squareSum_Y = 0;

    for (int i = 0; i < n; i++)
    {
        // sum of elements of array X.
        sum_X = sum_X + X[i];

        // sum of elements of array Y.
        sum_Y = sum_Y + Y[i];

        // sum of X[i] * Y[i].
        sum_XY = sum_XY + X[i] * Y[i];

        // sum of square of array elements.
        squareSum_X = squareSum_X +
                      X[i] * X[i];
        squareSum_Y = squareSum_Y +
                      Y[i] * Y[i];
    }

    // use formula for calculating
    // correlation coefficient.
    float corr = (float)(n * sum_XY -
                         sum_X * sum_Y) /
                 sqrt((n * squareSum_X -
                       sum_X * sum_X) *
                      (n * squareSum_Y -
                       sum_Y * sum_Y));

    return corr;
}

