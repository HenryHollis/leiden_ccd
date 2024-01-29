#include <igraph/igraph.h>
#include "GraphHelper.h"
#include "Optimiser.h"
#include "ModularityVertexPartition.h"
#include "ccdModularityVertexPartition.h"
#include <cstdio>

using std::cout;
using std::endl;

int main()
{
    igraph_t g;
    igraph_famous(&g, "Zachary");

    Graph graph(&g);

    ModularityVertexPartition part(&graph); //creates a ModularityVertexPartition called part
    Optimiser o; //create optimiser o
    std::cout<< "partition1 quality before optimization: " << part.quality() <<endl;
    o.optimise_partition(&part); //coptimise our ModularityVertexPartition obj
    cout << "Node\tCommunity" << endl;
    for (int i = 0; i < graph.vcount(); i++)
        cout << i << "\t" << part.membership(i) << endl;

    std::cout<< "quality AFTER optimization: " << part.quality() <<endl;

    igraph_destroy(&g);
}
