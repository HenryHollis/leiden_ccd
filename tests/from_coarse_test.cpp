//
// Created by Henry Hollis on 3/15/24.
//
#include <igraph/igraph.h>
#include "GraphHelper.h"
#include "Optimiser.h"
#include "LouvainOptimiser.h"
#include "ModularityVertexPartition.h"
#include "ccdModularityVertexPartition.h"
#include <cstdio>
#include <random>


using std::cout;
using std::endl;

vector<double> emat = { 0.41715202, 0.23533709, 0.7158693 , 0.76774877, 0.67170842,
                        0.85459821, 0.43317033, 0.76930908, 0.53764885, 0.42389998, 0.68937743,
                        0.20011722,
                        0.28176733, 0.24876795, 0.41850602, 0.9433223 , 0.65084243,
                        0.32621739,
                        0.09598426, 0.42357955, 0.72969605, 0.15087005, 0.6675312 ,
                        0.98522932,
                        0.8590361 , 0.20656894, 0.2966243 , 0.58040359, 0.66213785,
                        0.11772633,
                        0.02125281, 0.01135435, 0.19239381, 0.83530971, 0.68056444,
                        0.7165933 ,
                        0.15183117, 0.53466693, 0.54793279, 0.15587421, 0.38399057,
                        0.12140239,
                        0.39562934, 0.40733329, 0.00463059, 0.0048643 , 0.36128227,
                        0.90301726,
                        0.80607354, 0.33767792, 0.50406958, 0.83607252, 0.99133673,
                        0.86028985,
                        0.28966774, 0.13376067, 0.84430879, 0.3790256 , 0.88808115,
                        0.54172877,
                        0.32126778, 0.66190695, 0.56181687, 0.48214338, 0.36219704,
                        0.75276438,
                        0.95009022, 0.46247351, 0.77223523, 0.9479285 , 0.66446316,
                        0.4040721 };

vector<double> refmat = { 1.        , -0.92660504,  0.99832331, -0.59621204, -0.91019751,
                          0.47477778, -0.30982677,  0.75702281,  0.44129599, -0.00534183,
                          -0.84198001, -0.04100092,
                          -0.92660504,  1.        , -0.94681794,  0.25056121,  0.99913813,
                          -0.77088294, -0.07044538, -0.94715875, -0.0714668 ,  0.38098043,
                          0.98305764,  0.41371148,
                          0.99832331, -0.94681794,  1.        , -0.54874141, -0.93264553,
                          0.5249259 , -0.25427143,  0.79357436,  0.38861305, -0.06321621,
                          -0.87179728, -0.09876766,
                          -0.59621204,  0.25056121, -0.54874141,  1.        ,  0.21016028,
                          0.4235044 ,  0.94804477,  0.07321178, -0.98353211, -0.79963069,
                          0.06886645, -0.77770667,
                          -0.91019751,  0.99913813, -0.93264553,  0.21016028,  1.        ,
                          -0.79665886, -0.11179063, -0.95965709, -0.03000226,  0.41903067,
                          0.98981885,  0.45114512,
                          0.47477778, -0.77088294,  0.5249259 ,  0.4235044 , -0.79665886,
                          1.        ,  0.6896995 ,  0.93446851, -0.58025554, -0.88262934,
                          -0.87457811, -0.89883196,
                          -0.30982677, -0.07044538, -0.25427143,  0.94804477, -0.11179063,
                          0.6896995 ,  1.        ,  0.38669125, -0.98993046, -0.94912441,
                          -0.25209318, -0.93729033,
                          0.75702281, -0.94715875,  0.79357436,  0.07321178, -0.95965709,
                          0.93446851,  0.38669125,  1.        , -0.25225458, -0.65742302,
                          -0.98990684, -0.68387766,
                          0.44129599, -0.0714668 ,  0.38861305, -0.98353211, -0.03000226,
                          -0.58025554, -0.98993046, -0.25225458,  1.        ,  0.89499147,
                          0.11257201,  0.87851348,
                          -0.00534183,  0.38098043, -0.06321621, -0.79963069,  0.41903067,
                          -0.88262934, -0.94912441, -0.65742302,  0.89499147,  1.        ,
                          0.54399874,  0.99936387,
                          -0.84198001,  0.98305764, -0.87179728,  0.06886645,  0.98981885,
                          -0.87457811, -0.25209318, -0.98990684,  0.11257201,  0.54399874,
                          1.        ,  0.57357701,
                          -0.04100092,  0.41371148, -0.09876766, -0.77770667,  0.45114512,
                          -0.89883196, -0.93729033, -0.68387766,  0.87851348,  0.99936387,
                          0.57357701,  1.        };


int main(void) {
    igraph_t g;
    igraph_vector_int_t edges;

    /* Create a directed graph with no vertices or edges. */
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);

    /* Add 5 vertices. Vertex IDs will range from 0 to 4, inclusive. */
    igraph_add_vertices(&g, 6, NULL);

    /* Add 5 edges, specified as 5 consecutive pairs of vertex IDs
     * stored in an integer vector. */
    igraph_vector_int_init_int(&edges, 12,
                               0,1, 0,4, 4,1, 0,3, 3,2, 2,5);
    igraph_add_edges(&g, &edges, NULL);
    igraph_vector_int_destroy(&edges);
    igraph_write_graph_dot(&g, stdout);
    Graph graph(&g);
    ccdModularityVertexPartition part(&graph, {0,0,0,1,1,1}); //cre
//    part.setGeneSampleMatrix(emat, 12, 3);
    part.setGeneSampleMatrix(emat, 12, 6);

    part.setRefMatrix(refmat, 12, 12);

    cout << "Node\tCommunity" << endl;
    for (int i = 0; i < graph.vcount(); i++)
        cout << i << "\t" << part.membership(i) << endl;

    // Collapse graph (i.e. community graph)
    Graph* new_collapsed_graphs;
    Graph* current_graph = part.get_graph();
    vector<MutableVertexPartition*> current_partitions(1);
    current_partitions[0] = &part;
    MutableVertexPartition* new_collapsed_partition;

    new_collapsed_graphs = current_graph->collapse_graph(current_partitions[0]);
    // Create collapsed partition (i.e. default partition of each node in its own community).
    new_collapsed_partition = current_partitions[0]->create(new_collapsed_graphs);
    new_collapsed_partition->set_membership({2, 3});
    current_partitions[0]->from_coarse_partition(new_collapsed_partition);

    cout << "Node\tCommunity" << endl;
    for (int i = 0; i < graph.vcount(); i++)
        cout << i << "\t" << current_partitions[0]->membership(i) << endl;
    igraph_destroy(&g);

    return 0;
}