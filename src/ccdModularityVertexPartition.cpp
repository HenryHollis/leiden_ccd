#include "../include/ccdModularityVertexPartition.h"

#ifdef DEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif
ccdModularityVertexPartition::ccdModularityVertexPartition(Graph *graph, const vector<size_t> &membership,
                                                           const vector<std::vector<double>> &geneSampleMatrix):
                                                           MutableVertexPartition(graph, membership),
                                                           geneSampleMatrix(geneSampleMatrix){

}

ccdModularityVertexPartition::ccdModularityVertexPartition(Graph* graph,
                                                     vector<size_t> const& membership) :
        MutableVertexPartition(graph,
                               membership)
{ }

ccdModularityVertexPartition::ccdModularityVertexPartition(Graph* graph) :
        MutableVertexPartition(graph)
{ }

ccdModularityVertexPartition::~ccdModularityVertexPartition()
{ }

ccdModularityVertexPartition* ccdModularityVertexPartition::create(Graph* graph)
{
    return new ccdModularityVertexPartition(graph);
}

ccdModularityVertexPartition* ccdModularityVertexPartition::create(Graph* graph, vector<size_t> const& membership)
{
    return new ccdModularityVertexPartition(graph, membership);
}

ccdModularityVertexPartition *ccdModularityVertexPartition::create(Graph *graph, const vector<size_t> &membership,
                                                                   const vector<std::vector<double>> &geneSampleMatrix) {
    return new ccdModularityVertexPartition(graph, membership, geneSampleMatrix);
}

void ccdModularityVertexPartition::setGeneSampleMatrix(const vector<std::vector<double>> &geneSampleMatrix) {
    this->geneSampleMatrix = geneSampleMatrix;
}

const std::vector<std::vector<double>> &ccdModularityVertexPartition::getMatrix() {
    return geneSampleMatrix;
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
  double total_weight = this->graph->total_weight()*(2.0 - this->graph->is_directed());
  if (total_weight == 0.0)
    return 0.0;
  if (new_comm != old_comm)
  {
    #ifdef DEBUG
      cerr << "\t" << "old_comm: " << old_comm << endl;
    #endif
    double w_to_old = this->weight_to_comm(v, old_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_to_old: " << w_to_old << endl;
    #endif
    double w_from_old = this->weight_from_comm(v, old_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_from_old: " << w_from_old << endl;
    #endif
    double w_to_new = this->weight_to_comm(v, new_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_to_new: " << w_to_new << endl;
    #endif
    double w_from_new = this->weight_from_comm(v, new_comm);
    #ifdef DEBUG
      cerr << "\t" << "w_from_new: " << w_from_new << endl;
    #endif
    double k_out = this->graph->strength(v, IGRAPH_OUT);
    #ifdef DEBUG
      cerr << "\t" << "k_out: " << k_out << endl;
    #endif
    double k_in = this->graph->strength(v, IGRAPH_IN);
    #ifdef DEBUG
      cerr << "\t" << "k_in: " << k_in << endl;
    #endif
    double self_weight = this->graph->node_self_weight(v);
    #ifdef DEBUG
      cerr << "\t" << "self_weight: " << self_weight << endl;
    #endif
    double K_out_old = this->total_weight_from_comm(old_comm);
    #ifdef DEBUG
      cerr << "\t" << "K_out_old: " << K_out_old << endl;
    #endif
    double K_in_old = this->total_weight_to_comm(old_comm);
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
  return diff/m;
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

double ccdModularityVertexPartition::calcCCDsimple(const vector<std::vector<double>> &ref,
                                                   const vector<std::vector<double>> &emat, const string &method,
                                                   bool scale) {


    // Placeholder code (replace with actual implementation)
    size_t numRowsRef = ref.size();
    size_t numRowsEmat = emat.size();

    // Loop over matrix elements (simplified)
    for (size_t i = 0; i < numRowsRef; ++i) {
        for (size_t j = 0; j < numRowsEmat; ++j) {

        }
    }

    if (!ref.empty() && !ref[0].empty()) {
        // Number of columns is the size of any row (assuming all rows have the same size)
        size_t numColumns = ref[0].size();
        if (scale) {
            size_t nPairs = choose(numColumns, 2);
            ccd /= static_cast<double>(nPairs);
        }
    } else {
        std::cerr << "Matrix ref is empty or has empty rows." << std::endl;
    }

    // Placeholder: Scale ccd if requested


    return ccd;

}

// Function to calculate the factorial of a number

//calcDist = function(r1, r2) sqrt(sum((r1 - r2)^2, na.rm = TRUE))


long choose(size_t n, int k) {
    if (0 == k)
        return 1;
    return (n * choose(n - 1, k - 1)) / k;
}
/**
 * Spearman Correlation Code
 * Inspired from https://www.geeksforgeeks.org/program-spearmans-rank-correlation/
 */


// Utility Function to print vector
void printVector(const Vector &X) {
    for (auto i: X)
        std::cout <<i<<" ";
    std::cout << endl;
}

// Function returns the rank vector
// of the set of observations
Vector rankify(Vector &X) {

    int N = X.size();

    // Rank Vector
    Vector Rank_X(N);

    for(int i = 0; i < N; i++)
    {
        int r = 1, s = 1;

        // Count no of smaller elements
        // in 0 to i-1
        for(int j = 0; j < i; j++) {
            if (X[j] < X[i] ) r++;
            if (X[j] == X[i] ) s++;
        }

        // Count no of smaller elements
        // in i+1 to N-1
        for (int j = i+1; j < N; j++) {
            if (X[j] < X[i] ) r++;
            if (X[j] == X[i] ) s++;
        }

        // Use Fractional Rank formula
        // fractional_rank = r + (n-1)/2
        Rank_X[i] = r + (s-1) * 0.5;
    }

    // Return Rank Vector
    return Rank_X;
}

// function that returns
// Pearson correlation coefficient.
float cor(Vector &X, Vector &Y) {
    int n = X.size();
    float sum_X = 0, sum_Y = 0,
            sum_XY = 0;
    float squareSum_X = 0,
            squareSum_Y = 0;

    for (int i = 0; i < n; i++)
    {
        // sum of elements of array X.
        sum_X = sum_X + X[i];

        // sum of elements of array Y.
        sum_Y = sum_Y + Y[i];

        // sum of X[i] * Y[i].
        sum_XY = sum_XY + X[i] * Y[i];

        // sum of square of array elements.
        squareSum_X = squareSum_X +
                      X[i] * X[i];
        squareSum_Y = squareSum_Y +
                      Y[i] * Y[i];
    }

    // use formula for calculating
    // correlation coefficient.
    float corr = (float)(n * sum_XY -
                         sum_X * sum_Y) /
                 sqrt((n * squareSum_X -
                       sum_X * sum_X) *
                      (n * squareSum_Y -
                       sum_Y * sum_Y));

    return corr;
}








