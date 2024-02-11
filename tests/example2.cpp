#include <igraph/igraph.h>
#include "GraphHelper.h"
#include "Optimiser.h"
#include "ModularityVertexPartition.h"
#include "ccdModularityVertexPartition.h"
#include "LouvainOptimiser.h"
#include <cstdio>
#include <random>


using std::cout;
using std::endl;

int main() {
    // Set a fixed seed for reproducibility
    std::random_device rd;
    std::mt19937 gen(42);  // Use std::mt19937 with a random seed from std::random_device

    // Define the dimensions of the matrix
    const size_t rows = 12;
    const size_t cols = 34;
    igraph_t g;
    igraph_famous(&g, "Zachary");
    Graph graph(&g);

    // Initialize the matrix with zeros
    std::vector<std::vector<double>> geneSampleMatrix(rows, std::vector<double>(cols, 0.0));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            geneSampleMatrix[i][j] = std::generate_canonical<double, 10>(gen);// Example: random doubles between 0 and 1
        }
    }


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

    LouvainOptimiser o; //create optimiser o
    o.optimise_partition(&part); //coptimise our ModularityVertexPartition obj

    cout << "Node\tCommunity" << endl;
    for (int i = 0; i < graph.vcount(); i++)
        cout << i << "\t" << part.membership(i) << endl;

    std::cout<< "quality AFTER optimization: " << part.quality() <<endl;

    igraph_destroy(&g);
    return 0;
}