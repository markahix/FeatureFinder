#ifndef BASEDATATYPES_H
#define BASEDATATYPES_H
#include <vector>
#include <unordered_map>
#include <set>
#include <functional>


/*
     Type Definitions ONLY
*/

// Geometry
struct xyzCoords
{
    double x;
    double y;
    double z;
};

// Adjacency Matrix
typedef std::vector<std::vector<int>> AdjMat_int;


/*
    * INTERNAL HELPER TYPES
    * ---------------------
    * The public function signatures stay the same. Internally, however, we stop
    * solving subgraph matching by brute-force “choose k out of n” enumeration.
    *
    * Instead, we convert the legacy adjacency-matrix representation into a sparse
    * adjacency-list view and run a backtracking matcher with strong pruning.
    *
    * This is the same algorithmic idea that made the Python prototype fast:
    *   1. build candidate domains for each query atom
    *   2. choose the most constrained query atom next
    *   3. verify edge consistency incrementally
    *   4. only explore plausible mappings
    */

struct SearchGraph
{
    std::vector<int> labels;                          // encoded node labels: atomic_number * 10 + degree
    std::vector<std::vector<int>> neighbors;          // sparse adjacency list
    std::vector<std::unordered_map<int,int>> edge_label; // edge_label[u][v] = bond order label
    std::vector<int> degree;
};

struct PreparedQuery
{
    std::vector<int> orig_index_of_explicit;         // explicit (non-R) query vertices only
    std::vector<int> encoded_label;                  // encoded labels for explicit query atoms
    std::vector<std::vector<int>> neighbors;         // explicit query connectivity
    std::vector<std::unordered_map<int,int>> edge_label; // query edge_label[u][v] = bond order label
    std::vector<int> required_external_neighbors;    // count of attached R-group leaves for each explicit query atom
    std::vector<int> total_required_degree;          // exact degree required = explicit degree + wildcard attachments
};


#endif