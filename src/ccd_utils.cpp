//
// Created by Henry Hollis on 1/28/24.
//

#include "../include/ccd_utils.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
//Definition of the constant matrix
//const std::vector<double> ccd_utils::refCor = {
//        1.0000000,   0.77547090,  0.72492855,  0.27817942, -0.63637681, -0.60375141, -0.8614806, -0.7471112, -0.59455286, -0.8234182, -0.9146447, -0.8473980, //ARNTL
//        0.7754709,   1.00000000,  0.63439613,  0.07402797, -0.62632300, -0.34987550, -0.7461844, -0.6450780, -0.70865725, -0.7845410, -0.7654845, -0.7983427, //NPAS2
//        0.7249286,   0.63439613,  1.00000000,  0.06541974, -0.59727560, -0.30024636, -0.6031795, -0.6364953, -0.56958405, -0.7144612, -0.6455111, -0.7595101, //CLOCK
//        0.2781794,   0.07402797,  0.06541974,  1.00000000, -0.01245765, -0.72253596, -0.4099044, -0.1411756,  0.25538496, -0.0252816, -0.3401805, -0.0781101, //CRY1
//        -0.6363768, -0.62632300, -0.59727560, -0.01245765,  1.00000000,  0.28367324,  0.6234166,  0.6454257,  0.59510653,  0.6712806,  0.6618797,  0.7597038, //CRY2
//        -0.6037514, -0.34987550, -0.30024636, -0.72253596,  0.28367324,  1.00000000,  0.6772739,  0.4242223, -0.06776682,  0.3366267,  0.6955807,  0.3810191, //NR1D1
//        -0.8614806, -0.74618443, -0.60317949, -0.40990436,  0.62341661,  0.67727389,  1.0000000,  0.7132144,  0.52923596,  0.7673822,  0.9111478,  0.7487607, //NR1D2
//        -0.7471112, -0.64507795, -0.63649530, -0.14117556,  0.64542570,  0.42422234,  0.7132144,  1.0000000,  0.60794410,  0.7467579,  0.7732704,  0.7756198, //PER1
//        -0.5945529, -0.70865725, -0.56958405,  0.25538496,  0.59510653, -0.06776682,  0.5292360,  0.6079441,  1.00000000,  0.7868302,  0.5543211,  0.7530874, //PER2
//        -0.8234182, -0.78454102, -0.71446119, -0.02528160,  0.67128060,  0.33662668,  0.7673822,  0.7467579,  0.78683019,  1.0000000,  0.8117621,  0.8738338, //PER3
//        -0.9146447, -0.76548454, -0.64551113, -0.34018047,  0.66187971,  0.69558073,  0.9111478,  0.7732704,  0.55432112,  0.8117621,  1.0000000,  0.8443479, //DBP
//        -0.8473980, -0.79834269, -0.75951011, -0.07811010,  0.75970381,  0.38101906,  0.7487607,  0.7756198,  0.75308740,  0.8738338,  0.8443479,  1.0000000 //TEF
//};

double ccd_utils::calcCCDsimple(const std::vector<double> &ref, int num_ref_rows,
                     const std::vector<double> &emat, size_t emat_row, size_t emat_col,
                      bool scale) {

    std::vector<double> cormat = calcCorMat(emat, emat_row, emat_col); //calculate cormat of expression matrix
    if (num_ref_rows != emat_row || cormat.empty() || ref.empty() ) {
        throw std::invalid_argument("Matrices must be of the same size for calcCCDsimple");
    }
    double upperTriDiff = 0.0;
    //loop through triangular matrix and accumulate the difference between entries of ref and cormat
    for (size_t i = 0; i < num_ref_rows; ++i) {
        for (size_t j = i; j < num_ref_rows; ++j) {
            upperTriDiff += pow(ref[i * num_ref_rows + j] - cormat[i * num_ref_rows + j], 2);
        }
    }
    double ccd = sqrt(upperTriDiff);

    if (!ref.empty() ) {
        // Number of columns is the size of any row (assuming all rows have the same size)
        if (scale) {
            size_t nPairs = choose(num_ref_rows, 2);
            ccd /= static_cast<double>(nPairs);
        }
    } else
        std::cerr << "Matrix ref is empty or has empty rows." << std::endl;

    return ccd;
}

std::vector<double> ccd_utils::calcCorMat(const std::vector<double> &rect, size_t numRows, size_t numCols) {
   /* Takes rectangular matrix and calculated the gene x gene correlation matrix
    *
    */

    // Initialize the flat correlation matrix with zeros
    std::vector<double> correlationMatrix(numRows*numRows ,0.0);

    // Calculate pairwise correlations for each row
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = i; j < numRows; ++j) {
            if (i == j) {
                // Diagonal elements (correlation with itself) are always 1
                correlationMatrix[i * numRows + j] = 1.0;
            } else {
                // Pass a single row to the function
                std::vector<double> rowToProcessi(rect.begin() + i * numCols, rect.begin() + (i + 1) * numCols);
                std::vector<double> rowToProcessj(rect.begin() + j * numCols, rect.begin() + (j + 1) * numCols);
//                correlationMatrix[i * numRows + j] = correlationMatrix[i * numRows + j] = cor(rowToProcessi, rowToProcessj);

                std::vector<double> rank_x = rankVector(rowToProcessi);
                std::vector<double> rank_y = rankVector(rowToProcessj);
                correlationMatrix[i * numRows + j] = correlationMatrix[i * numRows + j] = cor(rank_x, rank_y);
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
// Function returns the rank vector
// of the set of observations
std::vector<double> ccd_utils::rankVector(const std::vector<double> &X) {

    int N = X.size();

    // Rank Vector
    std::vector<double> Rank_X(N);

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
double ccd_utils::cor(const std::vector<double> &X, const std::vector<double> &Y) {
    int n = X.size();
    double sum_X = 0, sum_Y = 0,
            sum_XY = 0;
    double squareSum_X = 0,
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
    double corr = (double)(n * sum_XY -
                         sum_X * sum_Y) /
                 sqrt((n * squareSum_X -
                       sum_X * sum_X) *
                      (n * squareSum_Y -
                       sum_Y * sum_Y));

    return corr;
}




std::vector<double> ccd_utils::sliceColumns(const std::vector<double> &matrix, const std::vector<size_t> &columnsToAccess,
                                            size_t nrow, size_t ncol)  {
    if (matrix.empty() || columnsToAccess.empty()) {
        std::cerr << "Matrix is empty or columns to slice are empty." << std::endl;
        return {}; // Return an empty matrix if either the matrix or columns are empty
    }

     // Create a vector to store the elements of the specified columns
    std::vector<double> columnsToProcess;
    // Copy the elements of the specified columns into the processing vector
//    for (int col : columnsToAccess) {
//        for (int row = 0; row < nrow; ++row) {
//            columnsToProcess.push_back(matrix[row * ncol + col]);
//        }
//    }

    // Copy the elements of the specified columns into the processing vector
    for (int row = 0; row < nrow; ++row) {
        for (int col : columnsToAccess) {
            columnsToProcess.push_back(matrix[row * ncol + col]);
        }
    }


    return columnsToProcess;
}
