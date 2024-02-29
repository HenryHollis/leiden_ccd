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
std::vector<double> refCor = {
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
int main() {
    // Set a fixed seed for reproducibility
    std::random_device rd;
    std::mt19937 gen(42);  // Use std::mt19937 with a random seed from std::random_device

    // Define the dimensions of the matrix
    const size_t rows = 12;
    const size_t cols = 100;
    igraph_t g;
    igraph_erdos_renyi_game_gnp(
            &g, 100, .01,
            false, false);
    Graph graph(&g);
    igraph_write_graph_dot(&g, stdout);
    // Initialize the matrix with zeros
    std::vector<double> geneSampleMatrix(rows*cols, 0.);

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            geneSampleMatrix[i * cols + j] = std::generate_canonical<double, 10>(gen);// Example: random doubles between 0 and 1
        }
    }


    //Set the partition matrix to the modified matrix
    ccdModularityVertexPartition part(&graph); //creates a ModularityVertexPartition called part
    part.setGeneSampleMatrix(geneSampleMatrix, 12, 100);
    part.setRefMatrix(refCor, 12, 12);
//    const std::vector<double>& matrix = part.getGeneMatrix();
//
//    // Print the modified matrix
//    std::cout << "\npart.matrix Matrix:" << std::endl;
//    for (size_t i = 0; i < rows; ++i) {
//        for (size_t j = 0; j < cols; ++j) {
//            std::cout << matrix[i * cols + j] << " ";
//        }
//        std::cout << std::endl;
//    }

    Optimiser o; //create optimiser o
    o.optimise_partition(&part); //optimise our ModularityVertexPartition obj

    cout << "Node\tCommunity" << endl;
    for (int i = 0; i < graph.vcount(); i++)
        cout << i << "\t" << part.membership(i) << endl;

    std::cout<< "quality AFTER optimization: " << part.quality() <<endl;

    igraph_destroy(&g);
    return 0;
}