#include <igraph/igraph.h>
#include "GraphHelper.h"
#include "Optimiser.h"
#include "ModularityVertexPartition.h"
#include "ccdModularityVertexPartition.h"
#include <cstdio>

using std::cout;
using std::endl;

int main() {
    // Define the dimensions of the matrix
    const size_t rows = 3;
    const size_t cols = 4;
    igraph_t g;
    igraph_famous(&g, "Zachary");
    Graph graph(&g);

    // Initialize the matrix with zeros
    std::vector<std::vector<double>> geneSampleMatrix(rows, std::vector<double>(cols, 0.0));

    // Modify some values in the matrix
    geneSampleMatrix[1][2] = 1.5;
    geneSampleMatrix[2][3] = 2.3;

    //Set the partition matrix to the modified matrix
    ccdModularityVertexPartition part(&graph); //creates a ModularityVertexPartition called part
    part.setGeneSampleMatrix(geneSampleMatrix);

    const std::vector<std::vector<double>>& matrix = part.getMatrix();

    // Print the modified matrix
    std::cout << "\npart.matrix Matrix:" << std::endl;
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}