//
// Created by Henry Hollis on 1/27/24.
//

#ifndef LEIDEN_CCD_CCDMODULARITYVERTEXPARTITION_H
#define LEIDEN_CCD_CCDMODULARITYVERTEXPARTITION_H
#define CCD_COMM_SIZE 2
#include <unordered_map>
#include <numeric>
#include <functional>


#include "MutableVertexPartition.h"
class LIBLEIDENALG_EXPORT ccdModularityVertexPartition : public MutableVertexPartition
{
public:
    //constructors
    // Constructor with matrix initialization
    ccdModularityVertexPartition(Graph* graph,
                                 vector<size_t> const& membership,
                                 const std::vector<double>& geneSampleMatrix);
    ccdModularityVertexPartition(Graph* graph,
                                 vector<size_t> const& membership);
    ccdModularityVertexPartition(Graph* graph);

    virtual ~ccdModularityVertexPartition();
    virtual ccdModularityVertexPartition* create(Graph* graph);
    virtual ccdModularityVertexPartition* create(Graph* graph, vector<size_t> const& membership);
    virtual ccdModularityVertexPartition* create(Graph* graph, vector<size_t> const& membership, const std::vector<double>& geneSampleMatrix);

    virtual double diff_move(size_t v, size_t new_comm);
    virtual double quality();
    //  Setter method for the matrix
    void setGeneSampleMatrix(const std::vector<double>& geneSampleMatrix, size_t rows, size_t cols);
    void setRefMatrix(const std::vector<double>& refMat, size_t rows, size_t cols);

    // Getter for geneSampleMatrix
    const std::vector<double>& getGeneMatrix();
    const std::vector<double> & getRefMatrix();

    void set_fine_membership(vector<size_t> const& new_membership);
    void from_coarse_partition_fine(vector<size_t> const& coarse_partition_membership, vector<size_t> const& coarse_node);
    void relabel_communities(vector<size_t> const& new_comm_id) override;
    void move_node(size_t v,size_t new_comm) override;
protected:
    vector<size_t> _fine_membership; // Membership vector, that never is collapse
    std::unordered_map<size_t, vector<size_t>> _membership_dict;
    static static std::unordered_map<size_t, std::vector<size_t>> membership_to_dict(vector<size_t> &membership);

private:
    // Matrix representing genes and samples
    std::vector<double> geneSampleMatrix;
    size_t geneMatRows;
    size_t geneMatCols;
    std::vector<double> refMatrix;
    size_t refMatRows;
    size_t refMatCols;
    struct vecHash {
        size_t operator()(const std::vector<size_t>& v) const;
    };
    std::unordered_map<std::vector<size_t>, double, vecHash> ccdCache;
};

#endif //LEIDEN_CCD_CCDMODULARITYVERTEXPARTITION_H
