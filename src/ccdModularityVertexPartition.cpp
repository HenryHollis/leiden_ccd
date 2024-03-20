#include "ccdModularityVertexPartition.h"
#include "ccd_utils.h"

#ifdef DEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif
ccdModularityVertexPartition::ccdModularityVertexPartition(Graph *graph, const vector<size_t> &membership,
                                                           const vector<double> &geneSampleMatrix):
        MutableVertexPartition(graph, membership),
        geneSampleMatrix(geneSampleMatrix){
}

ccdModularityVertexPartition::ccdModularityVertexPartition(Graph* graph,
                                                           vector<size_t> const& membership) :
        MutableVertexPartition(graph, membership)
{}

ccdModularityVertexPartition::ccdModularityVertexPartition(Graph* graph) :
        MutableVertexPartition(graph)
{}

ccdModularityVertexPartition::~ccdModularityVertexPartition()
{ }

ccdModularityVertexPartition* ccdModularityVertexPartition::create(Graph* graph)
{
    std::vector<double> GeneMatrix = this->geneSampleMatrix;
    size_t geneMatRows = this-> geneMatRows;
    size_t geneMatCols = this-> geneMatCols;

    std::vector<double> RefMat = this->refMatrix;
    size_t refMatRows = this-> refMatRows;
    size_t refMatCols = this-> refMatCols;

    auto* tmp = new ccdModularityVertexPartition(graph);
    tmp->geneSampleMatrix = GeneMatrix;
    tmp->refMatrix = RefMat;
    tmp->refMatRows = refMatRows;
    tmp->refMatCols = refMatCols;
    tmp->geneMatRows = geneMatRows;
    tmp->geneMatCols = geneMatCols;
    tmp->_fine_membership = this->_fine_membership;
    tmp->_membership_dict = membership_to_dict(this->_fine_membership);
    return tmp;

}

ccdModularityVertexPartition* ccdModularityVertexPartition::create(Graph* graph, vector<size_t> const& membership)
{
    std::vector<double> GeneMatrix = this->geneSampleMatrix;
    size_t geneMatRows = this-> geneMatRows;
    size_t geneMatCols = this-> geneMatCols;

    std::vector<double> refMat = this->refMatrix;
    size_t refMatRows = this-> refMatRows;
    size_t refMatCols = this-> refMatCols;

    auto* tmp = new  ccdModularityVertexPartition(graph, membership);
    tmp->geneSampleMatrix = GeneMatrix;
    tmp->refMatrix = refMat;
    tmp->refMatRows = refMatRows;
    tmp->refMatCols = refMatCols;
    tmp->geneMatRows = geneMatRows;
    tmp->geneMatCols = geneMatCols;

    return tmp;
}

ccdModularityVertexPartition *ccdModularityVertexPartition::create(Graph *graph, const vector<size_t> &membership,
                                                                   const vector<double> &geneSampleMatrix) {
    return new ccdModularityVertexPartition(graph, membership, geneSampleMatrix);
}

void ccdModularityVertexPartition::setGeneSampleMatrix(const vector<double> &geneSampleMatrix, size_t rows, size_t cols) {
    if( !geneSampleMatrix.empty()){
        this->geneSampleMatrix = geneSampleMatrix;
        this->geneMatRows = rows;
        this->geneMatCols = cols;
    } else
        throw std::invalid_argument("Gene expression must not be empty");
}

void ccdModularityVertexPartition::setRefMatrix(const vector<double> &refMat, size_t rows, size_t cols) {
    if( !refMat.empty()){
        this->refMatrix = refMat;
        this->refMatRows = rows;
        this->refMatCols = cols;
    } else
        throw std::invalid_argument("Reference Matrix must not be empty");
}

const std::vector<double> &ccdModularityVertexPartition::getGeneMatrix() {
    return geneSampleMatrix;
}
const std::vector<double> &ccdModularityVertexPartition::getRefMatrix() {
    return refMatrix;
}


size_t ccdModularityVertexPartition::vecHash::operator()(const std::vector<size_t>& v) const {
    return std::accumulate(v.begin(), v.end(), 0, [](size_t a, size_t b) { return a ^ (b + 0x9e3779b9 + (a << 6) + (a >> 2)); });
}

std::unordered_map<size_t, std::vector<size_t>> ccdModularityVertexPartition::membership_to_dict(vector<size_t> &membership) {

    std::unordered_map<size_t, std::vector<size_t>> myMap;
    // Inserting unique values from A into myMap along with their indices
    for (size_t i = 0; i < membership.size(); ++i) {
        myMap[membership[i]].push_back(i);
    }

    return myMap;

}

void ccdModularityVertexPartition::set_fine_membership(const vector<size_t> &new_membership) {
  this->_fine_membership = new_membership;
  this->_membership_dict = membership_to_dict(this->_membership);
}

void ccdModularityVertexPartition::move_node(size_t v, size_t new_comm) {

    size_t old_comm = this->_membership[v]; //what community is v in?
    std::vector<size_t> nodes_in_v; // v itself has nodes (after partition collapse).
    if (this->_membership_dict.find(old_comm) != this->_membership_dict.end()) {
        nodes_in_v = this->_membership_dict.find(old_comm)->second;
    }else{
        std::cerr << "Key " << old_comm << " not found in _membership_dict." << std::endl;
    }

    // Iterate over the indices of nodes_in_v and update corresponding elements in _fine_membership
    for (int index : nodes_in_v) {
        if (index >= 0 && index < this->_fine_membership.size()) {
            this->_fine_membership[index] = new_comm;
        } else
            std::cerr << "Error: Index out of bounds." << std::endl;
    }

    this->_membership_dict = membership_to_dict(this->_fine_membership);
    MutableVertexPartition::move_node(v, new_comm);

}

void ccdModularityVertexPartition::relabel_communities(const vector<size_t> &new_comm_id){
    if (this->_n_communities != new_comm_id.size()) {
        throw Exception("Problem swapping community labels. Mismatch between n_communities and new_comm_id vector.");
    }

    size_t n = this->graph->vcount();

    for (size_t i = 0; i < n; i++){
        this->_membership[i] = new_comm_id[this->_membership[i]];
        //this->_fine_membership[i] = this->_membership[i];
    }
    this->from_coarse_partition_fine(this->membership(), this->_fine_membership);
    this->_membership_dict = membership_to_dict(this->_fine_membership);
    this->update_n_communities();
    size_t nbcomms = this->n_communities();

    vector<double> new_total_weight_in_comm(nbcomms, 0.0);
    vector<double> new_total_weight_from_comm(nbcomms, 0.0);
    vector<double> new_total_weight_to_comm(nbcomms, 0.0);
    vector<double> new_csize(nbcomms, 0);
    vector<size_t> new_cnodes(nbcomms, 0);

    // Relabel community admin
    for (size_t c = 0; c < new_comm_id.size(); c++) {
        size_t new_c = new_comm_id[c];
        if (this->_cnodes[c] > 0) {
            new_total_weight_in_comm[new_c] = this->_total_weight_in_comm[c];
            new_total_weight_from_comm[new_c] = this->_total_weight_from_comm[c];
            new_total_weight_to_comm[new_c] = this->_total_weight_to_comm[c];
            new_csize[new_c] = this->_csize[c];
            new_cnodes[new_c] = this->_cnodes[c];
        }
    }

    this->_total_weight_in_comm = new_total_weight_in_comm;
    this->_total_weight_from_comm = new_total_weight_from_comm;
    this->_total_weight_to_comm = new_total_weight_to_comm;
    this->_csize = new_csize;
    this->_cnodes = new_cnodes;

    this->_empty_communities.clear();
    for (size_t c = 0; c < nbcomms; c++) {
        if (this->_cnodes[c] == 0) {
            this->_empty_communities.push_back(c);
        }
    }

    // invalidate cached weight vectors
    for (size_t c : this->_cached_neigh_comms_from)
        this->_cached_weight_from_community[c] = 0;
    this->_cached_neigh_comms_from.clear();
    this->_cached_weight_from_community.resize(nbcomms, 0);
    this->_current_node_cache_community_from = n + 1;

    for (size_t c : this->_cached_neigh_comms_to)
        this->_cached_weight_to_community[c] = 0;
    this->_cached_neigh_comms_to.clear();
    this->_cached_weight_to_community.resize(nbcomms, 0);
    this->_current_node_cache_community_to = n + 1;

    for (size_t c : this->_cached_neigh_comms_all)
        this->_cached_weight_all_community[c] = 0;
    this->_cached_neigh_comms_all.clear();
    this->_cached_weight_all_community.resize(nbcomms, 0);
    this->_current_node_cache_community_all = n + 1;

#ifdef DEBUG
    if (this->_csize.size() < this->_n_communities ||
        this->_cnodes.size() < this->_n_communities ||
        this->_total_weight_in_comm.size() < this->_n_communities ||
        this->_total_weight_to_comm.size() < this->_n_communities ||
        this->_total_weight_from_comm.size() < this->_n_communities ||
        this->_cached_weight_from_community.size() < this->_n_communities ||
        this->_cached_weight_to_community.size() < this->_n_communities ||
        this->_cached_weight_all_community.size() < this->_n_communities) {
      cerr << "ERROR: MutableVertexPartition bookkeeping is too small after rearrange_community_labels." << endl;
    }

    this->init_admin();

    for (size_t c = 0; c < this->_n_communities; c++) {
      if (fabs(new_total_weight_in_comm[c] - this->_total_weight_in_comm[c]) > 1e-6 ||
          fabs(new_total_weight_from_comm[c] - this->_total_weight_from_comm[c]) > 1e-6 ||
          fabs(new_total_weight_to_comm[c] - this->_total_weight_to_comm[c]) > 1e-6 ||
          new_csize[c] != this->_csize[c] ||
          new_cnodes[c] != this->_cnodes[c]) {
        cerr << "ERROR: MutableVertexPartition bookkeeping is incorrect after rearrange_community_labels." << endl;
        cerr << "Community c has " << endl
             << "total_weight_in_comm=" << new_total_weight_in_comm[c]
             << " (should be " << this->_total_weight_in_comm[c] << ")" << endl
             << "total_weight_from_comm=" << new_total_weight_from_comm[c]
             << " (should be " << this->_total_weight_from_comm[c] << ")" << endl
             << "total_weight_to_comm=" << new_total_weight_to_comm[c]
             << " (should be " << this->_total_weight_to_comm[c] << ")" << endl
             << "csize=" << new_csize[c]
             << " (should be " << this->_csize[c] << ")" << endl
             << "cnodes=" << new_cnodes[c]
             << " (should be " << this->_cnodes[c] << ")" << endl;
      }
    }
#endif
}


void ccdModularityVertexPartition::from_coarse_partition_fine(const vector<size_t> &coarse_partition_membership,
                                                         const vector<size_t> &coarse_node) {
    /*This does the same as from_coarse_partition in MutableVertexPartition class
     * except that it modifies the _fine_membership member not the membership member.
     * In short after modification of a collapsed graph, this ensures that _fine_membership updates.
     * example:
     *
     * Finest partition before                                  Finest Partition After
     * nodes: 0, 1, 2, 3     -----from_coarse_partition--->    nodes: 0, 1, 2, 3
     * comm:  1, 0, 2, 1                                       comm:  1, 0, 0, 1
     *
     * Coarser Partition:
     * node = 0, 1, 2
     * comm = 0, 1, 0
     */

    // Read the coarser partition
    for (size_t v = 0; v < coarse_node.size(); v++)
    {
        // In the coarser partition, the node should have the community id
        // as represented by the coarser_membership vector
        size_t v_level2 = coarse_node[v];

        // In the coarser partition, this node is represented by v_level2
        size_t v_comm_level2 = coarse_partition_membership[v_level2];

        // Set local membership to community found for node at second level
        this->_fine_membership[v] = v_comm_level2;
    }

//    this->clean_mem();
//    this->init_admin();
}

/*****************************************************************************
  Returns the difference in modularity if we move a node to a new community
*****************************************************************************/
double ccdModularityVertexPartition::diff_move(size_t v, size_t new_comm)
{
#ifdef DEBUG
    cerr << "double ccdModularityVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
#endif
    size_t old_comm = this->_membership[v]; //what community is v in?
    double diff = 0.0;
    double ccd_diff = NAN;
    double total_weight = this->graph->total_weight()*(2.0 - this->graph->is_directed());

    std::vector<size_t> nodes_in_v; //NOTE: this could be different from nodes_in_old_comm_v when
                                    // v itself has nodes (after parition collapse).
    if (this->_membership_dict.find(old_comm) != this->_membership_dict.end()) {
        nodes_in_v = this->_membership_dict.find(old_comm)->second;
    }else{
        std::cerr << "Key " << old_comm << " not found in _membership_dict." << std::endl;
    }
    std::vector<size_t>  Nodes_in_old_comm_v;
    if (this->_membership_dict.find(old_comm) != this->_membership_dict.end()) {
        Nodes_in_old_comm_v = this->_membership_dict.find(old_comm)->second;
    }else{
            std::cerr << "Key " << old_comm << " not found in _membership_dict." << std::endl;
    }
    std::vector<size_t> Nodes_in_new_comm_no_v;
    if (this->_membership_dict.find(new_comm) != this->_membership_dict.end()) {
        Nodes_in_new_comm_no_v = this->_membership_dict.find(new_comm)->second;
    }
    std::vector<size_t> Nodes_in_new_comm_v;
    if (this->_membership_dict.find(new_comm) != this->_membership_dict.end()) {
        Nodes_in_new_comm_v = this->_membership_dict.find(new_comm)->second;
    }
    //Nodes_in_new_comm_v.push_back(v); // the nodes in the new community union v
    Nodes_in_new_comm_v.insert(Nodes_in_new_comm_v.end(), nodes_in_v.begin(), nodes_in_v.end());

    std::vector<size_t> Nodes_in_old_comm_no_v;
    if (this->_membership_dict.find(old_comm) != this->_membership_dict.end()) {
        Nodes_in_old_comm_no_v = this->_membership_dict.find(old_comm)->second;
    }else{
        std::cerr << "Key " << old_comm << " not found in _membership_dict." << std::endl;
    }

    // Define a lambda function to check if an element is in the array_to_delete
    auto is_in_array_to_delete = [&](int val) {
        return std::find(std::begin(nodes_in_v), std::end(nodes_in_v), val) != std::end(nodes_in_v);
    };

    // Use std::remove_if with the lambda function to remove elements from vec
    Nodes_in_old_comm_no_v.erase(std::remove_if(Nodes_in_old_comm_no_v.begin(), Nodes_in_old_comm_no_v.end(), is_in_array_to_delete), Nodes_in_old_comm_no_v.end());

    //Change in ccd should be [ccd(new+v) + ccd(old - v)] - [ccd(old + v) + ccd(new - v)]
    double old_ccd_v = 100;
    double new_ccd_no_v = 100;
    double old_ccd_no_v = 100;
    double new_ccd_w_v = 100;
    if (total_weight == 0.0)
        return 0.0;
    if (new_comm != old_comm)
    {
        // ********CALC CCD*************
        std::vector<double> emat = this->getGeneMatrix(); //Get the expression matrix associated with the partition object
        std::vector<double> refmat = this->getRefMatrix();
        // calculate ccd in old community if enough nodes are aggregated into c's community:
        if (CCD_COMM_SIZE < Nodes_in_old_comm_v.size()) {
            auto it = this->ccdCache.find(Nodes_in_old_comm_v);
            if (it != this->ccdCache.end()) {
                // Result is already in the cache, return it
                //        old_ccd =  it->second;
                //        std::cout << "old ccd found" <<std::endl;
            } else{
                //calculate the result and store it
                // try{
                //     std::vector<double> comm_emat = ccd_utils::sliceColumns(emat, Nodes_in_old_comm, this->geneMatRows, Nodes_in_old_comm.size() );
                std::vector<double> comm_emat_old = ccd_utils::sliceColumns(emat, Nodes_in_old_comm_v, this->geneMatRows, this->geneMatCols);

                old_ccd_v = ccd_utils::calcCCDsimple(refmat, this->refMatRows, comm_emat_old, this->geneMatRows,Nodes_in_old_comm_v.size(), false);
                //      this->ccdCache[Nodes_in_old_comm] = old_ccd;
                // }catch (const std::out_of_range& e) {
                //    std::cerr << "Exception caught: " << e.what() << std::endl;
                // }

            }
        }
        if (CCD_COMM_SIZE < Nodes_in_old_comm_no_v.size()) {
            std::vector<double> comm_emat_old2 = ccd_utils::sliceColumns(emat, Nodes_in_old_comm_no_v, this->geneMatRows, this->geneMatCols);
            old_ccd_no_v = ccd_utils::calcCCDsimple(refmat, this->refMatRows, comm_emat_old2, this->geneMatRows,Nodes_in_old_comm_no_v.size(), false);

        }
        if (CCD_COMM_SIZE < Nodes_in_new_comm_v.size()) {
            std::vector<double> comm_emat_new2 = ccd_utils::sliceColumns(emat,  Nodes_in_new_comm_v, this->geneMatRows, this->geneMatCols);
            new_ccd_w_v = ccd_utils::calcCCDsimple(refmat, this->refMatRows, comm_emat_new2, this->geneMatRows, Nodes_in_new_comm_v.size(), false);

        }
            //calc ccd of adding v into new community
        if (CCD_COMM_SIZE < Nodes_in_new_comm_no_v.size()) {
            auto it = this->ccdCache.find(Nodes_in_new_comm_no_v);
            if (it != this->ccdCache.end()) {
                // Result is already in the cache, return it
                //std::cout << "new ccd found" <<std::endl;
                // new_ccd = it->second;
            }else{
                //    calculate the result and store it
                //      try{
                std::vector<double> comm_emat_new = ccd_utils::sliceColumns(emat,  Nodes_in_new_comm_no_v, this->geneMatRows, this->geneMatCols);
                new_ccd_no_v = ccd_utils::calcCCDsimple(refmat, this->refMatRows, comm_emat_new, this->geneMatRows, Nodes_in_new_comm_no_v.size(), false);
                //    this->ccdCache[Nodes_in_new_comm] = new_ccd;
                //      }catch (const std::out_of_range& e) {
                //    std::cerr << "Exception caught: " << e.what() << std::endl;
                //      }


            }
        }
        //****************************
#ifdef DEBUG
        cerr << "\t" << "old_comm: " << old_comm << endl;
#endif
        double w_to_old = this->weight_to_comm(v, old_comm); //sum of edges going to old community
#ifdef DEBUG
        cerr << "\t" << "w_to_old: " << w_to_old << endl;
#endif
        double w_from_old = this->weight_from_comm(v, old_comm); //sum of edges coming from old community
#ifdef DEBUG
        cerr << "\t" << "w_from_old: " << w_from_old << endl;
#endif
        double w_to_new = this->weight_to_comm(v, new_comm); //sum of edges going to new community
#ifdef DEBUG
        cerr << "\t" << "w_to_new: " << w_to_new << endl;
#endif
        double w_from_new = this->weight_from_comm(v, new_comm); //sum of edges coming from new community
#ifdef DEBUG
        cerr << "\t" << "w_from_new: " << w_from_new << endl;
#endif
        double k_out = this->graph->strength(v, IGRAPH_OUT); //sum of all edges leaving node v
#ifdef DEBUG
        cerr << "\t" << "k_out: " << k_out << endl;
#endif
        double k_in = this->graph->strength(v, IGRAPH_IN); //sum of all edges coming into v
#ifdef DEBUG
        cerr << "\t" << "k_in: " << k_in << endl;
#endif
        double self_weight = this->graph->node_self_weight(v);
#ifdef DEBUG
        cerr << "\t" << "self_weight: " << self_weight << endl;
#endif
        double K_out_old = this->total_weight_from_comm(old_comm); //total weights of edges leaving OLD community
#ifdef DEBUG
        cerr << "\t" << "K_out_old: " << K_out_old << endl;
#endif
        double K_in_old = this->total_weight_to_comm(old_comm);  //total weights of edges Entering OLD community
#ifdef DEBUG
        cerr << "\t" << "K_in_old: " << K_in_old << endl;
#endif
        double K_out_new = this->total_weight_from_comm(new_comm) + k_out;
#ifdef DEBUG
        cerr << "\t" << "K_out_new: " << K_out_new << endl;
#endif
        double K_in_new = this->total_weight_to_comm(new_comm) + k_in;
#ifdef DEBUG
        cerr << "\t" << "K_in_new: " << K_in_new << endl;
      cerr << "\t" << "total_weight: " << total_weight << endl;
#endif
        double diff_old = (w_to_old - k_out*K_in_old/total_weight) + \
               (w_from_old - k_in*K_out_old/total_weight);
#ifdef DEBUG
        cerr << "\t" << "diff_old: " << diff_old << endl;
#endif
        double diff_new = (w_to_new + self_weight - k_out*K_in_new/total_weight) + \
               (w_from_new + self_weight - k_in*K_out_new/total_weight);
#ifdef DEBUG
        cerr << "\t" << "diff_new: " << diff_new << endl;
#endif
        diff = diff_new - diff_old;
#ifdef DEBUG
        cerr << "\t" << "diff: " << diff << endl;
#endif
        ccd_diff = (old_ccd_v + new_ccd_no_v) - (new_ccd_w_v + old_ccd_no_v); //negative number returns smaller score
        //ccd_diff = (old_ccd_v) - new_ccd_w_v;
    }
#ifdef DEBUG
    cerr << "exit double ccdModularityVertexPartition::diff_move((" << v << ", " << new_comm << ")" << endl;
    cerr << "return " << diff << endl << endl;
#endif
    double m;
    if (this->graph->is_directed())
        m = this->graph->total_weight();
    else
        m = 2*this->graph->total_weight();
    //ccd_diff = isfinite(ccd_diff) ? ccd_diff : 0.0;
    // Check if the result is within a tolerance of zero
//    const double tolerance = 1e-15;  // Adjust this threshold based on your needs

    double result = diff/m  + 0.1 * ccd_diff;
    std::cout <<"v: " << v<< "; new comm: " << new_comm <<"; old_com:" << old_comm <<"; old ccd w v:" << old_ccd_v <<"; old ccd no v:" << old_ccd_no_v  <<"; new_ccd_w_v:" <<  new_ccd_w_v << "; new_ccd_no_v:" << new_ccd_no_v<< "; ccd diff: "<<ccd_diff <<"; mod: "<<(diff/m) << "; res:" <<result << endl;

//if (std::abs(result) < tolerance) {
//        return 0.0;
//    }
    return result;
}


/*****************************************************************************
  Give the modularity of the partition.

  We here use the unscaled version of modularity, in other words, we don"t
  normalise by the number of edges.
******************************************************************************/
double ccdModularityVertexPartition::quality()
{
#ifdef DEBUG
    cerr << "double ccdModularityVertexPartition::quality()" << endl;
#endif
    double mod = 0.0;

    double m;
    if (this->graph->is_directed())
        m = this->graph->total_weight();
    else
        m = 2*this->graph->total_weight();

    if (m == 0)
        return 0.0;

    for (size_t c = 0; c < this->n_communities(); c++)
    {
        double w = this->total_weight_in_comm(c);
        double w_out = this->total_weight_from_comm(c);
        double w_in = this->total_weight_to_comm(c);
#ifdef DEBUG
        double csize = this->csize(c);
      cerr << "\t" << "Comm: " << c << ", size=" << csize << ", w=" << w << ", w_out=" << w_out << ", w_in=" << w_in << "." << endl;
#endif
        mod += w - w_out*w_in/((this->graph->is_directed() ? 1.0 : 4.0)*this->graph->total_weight());
    }
    double q = (2.0 - this->graph->is_directed())*mod;
#ifdef DEBUG
    cerr << "exit double ccdModularityVertexPartition::quality()" << endl;
    cerr << "return " << q/m << endl << endl;
#endif
    return q/m;
}













