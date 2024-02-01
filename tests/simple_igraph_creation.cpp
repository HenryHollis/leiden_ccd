#include <igraph/igraph.h>
#include "GraphHelper.h"
#include "Optimiser.h"
#include "ModularityVertexPartition.h"
#include "ccdModularityVertexPartition.h"
#include <cstdio>
#include <random>


using std::cout;
using std::endl;

int main(void) {
    igraph_t g;
    igraph_vector_int_t edges;

    /* Create a directed graph with no vertices or edges. */
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);

    /* Add 5 vertices. Vertex IDs will range from 0 to 4, inclusive. */
    igraph_add_vertices(&g, 7, NULL);

    /* Add 5 edges, specified as 5 consecutive pairs of vertex IDs
     * stored in an integer vector. */
    igraph_vector_int_init_int(&edges, 16,
                               0,1, 0,2, 2,1, 1,3, 3,4, 4,5,  5,3,  3,6);
    igraph_add_edges(&g, &edges, NULL);
    igraph_vector_int_destroy(&edges);
    igraph_write_graph_dot(&g, stdout);
    Graph graph(&g);
    ModularityVertexPartition part(&graph); //cre

    Optimiser o; //create optimiser o
    o.optimise_partition(&part); //coptimise our ModularityVertexPartition obj

    cout << "Node\tCommunity" << endl;
    for (int i = 0; i < graph.vcount(); i++)
        cout << i << "\t" << part.membership(i) << endl;


    igraph_destroy(&g);

    return 0;
}