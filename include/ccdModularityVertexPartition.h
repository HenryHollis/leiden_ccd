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

protected:
    vector<size_t> _fine_membership; // Membership vector, that never is collapse

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
