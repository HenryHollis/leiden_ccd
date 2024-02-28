#include <cmath>
#include <cstdio>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include "ccd_utils.h"
#include <chrono>
#include <random>


using namespace std;

int getHighestVal(int n, vector<double> arr) {
    int highest = 0;
    for (int i = 0; i < n; i++) {
        if (arr[i] > arr[highest])
            highest = i;
    }

    return highest;
}

vector<int> getRank(int n, vector<double> arr) {
    vector<int> rank(n);
    vector<bool> used(n);
    for (int i = 0; i < n; i++)
        used[i] = false;
    int lowestVal = getHighestVal(n, arr);

    for (int i = 1; i <= n; i++) {
        for (int j = 0; j < n; j++) {
            if (used[j] == false && arr[lowestVal] > arr[j])
                lowestVal = j;
        }

        rank[lowestVal] = i;
        used[lowestVal] = true;
        lowestVal = getHighestVal(n, arr);
    }

    return rank;
}

double diSquared(int n, vector<int> xRank, vector<int> yRank) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += pow(xRank[i] - yRank[i], 2);
    }

    return sum;
}

double spearman(int n, vector<int> xRank, vector<int> yRank) {
    return 1 - ((6 * diSquared(n, xRank, yRank)) / (n * (pow(n, 2) - 1)));
}

// Function to generate a vector of n random double values
std::vector<double> generateRandomVector(int n, double minVal, double maxVal) {
    // Seed for the random number generator
    std::random_device rd;

    // Mersenne Twister random number generator
    std::mt19937 gen(rd());

    // Uniform distribution between minVal and maxVal
    std::uniform_real_distribution<double> dis(minVal, maxVal);

    // Generate random values and fill the vector
    std::vector<double> randomVector;
    randomVector.reserve(n);
    for (int i = 0; i < n; ++i) {
        randomVector.push_back(dis(gen));
    }

    return randomVector;
}

int main() {
    // Set the size of the vector and the range of random values
    int n = 10000;
    double minValue = 0.0;
    double maxValue = 1.0;

    // Generate a vector of n random double values
    std::vector<double> x = generateRandomVector(n, minValue, maxValue);
    std::vector<double> y = generateRandomVector(n, minValue, maxValue);

    // Recording the timestamp at the start of the code
    auto beg = chrono::high_resolution_clock::now();
    cout << fixed << setprecision(3) << spearman(n, getRank(n, x), getRank(n, y)) << endl;
    // Taking a timestamp after the code is ran
    auto end = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::microseconds>(end - beg);
    // Displaying the elapsed time
    std::cout << "Elapsed Time: " << duration.count() <<"\n";
    auto beg1 = chrono::high_resolution_clock::now();

    std::vector<double> rank_x = ccd_utils::rankVector(x);
    std::vector<double> rank_y = ccd_utils::rankVector(y);
    cout << fixed << setprecision(3) << ccd_utils::cor(rank_x, rank_y) << endl;
    auto end1 = chrono::high_resolution_clock::now();
    auto duration1 = duration_cast<chrono::microseconds>(end1 - beg1);
    // Displaying the elapsed time
    std::cout << "Elapsed Time: " << duration1.count();
    return 0;
}