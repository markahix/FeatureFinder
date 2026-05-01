#ifndef SINGLEMOLECULE_H
#define SINGLEMOLECULE_H

#include "basedatatypes.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <set>
#include <ctime>
#include <cctype>

// DataType Definitions.

// Utilities
std::vector<std::string> ParsePDBLinesInOrder(std::string pdbfile); // reorder PDB in descending atomic mass order. (keeps matrices neater.)
void trim(std::string &s);
bool val_in_array(int val, std::vector<int> my_array);
bool Does_Vector_Contain_Vector(std::vector<int> bigger, std::vector<int> smaller);
bool any_val_in_array(std::set<int> vals, std::set<int>my_array);


// Geometry
const double RADIANS_TO_DEGREES = 57.2958;
double norm_of_vector(xyzCoords vector);
double calc_distance_between_points(xyzCoords pointA, xyzCoords pointB);


// Chemistry
int get_atomic_number(std::string symbol);
double get_covalent_radius(std::string symbol);
double get_covalent_radius(int number);
double get_vdw_radius(std::string symbol);
double get_vdw_radius(int number);
double get_atomic_mass(std::string symbol);
double get_atomic_mass(int number);

// Adjacency Matrices
class AdjacencyMatrix
{
    public:
        AdjacencyMatrix(std::vector<int> incoming_nodes, std::vector<std::pair<int,int>> incoming_edges);
        ~AdjacencyMatrix();
        std::string PrintAdjacencyMatrix();

    private:
        std::vector<int> node_weights;
        std::vector<std::pair<int,int>> edges;
        AdjMat_int adjacency_matrix;
};
AdjMat_int Generate_Adjacency_Matrix(std::vector<int> &atomic_numbers, std::vector<std::pair<int,int>> edges);
AdjMat_int Generate_Adjacency_Matrix(std::vector<int> &atomic_numbers, std::vector<std::pair<int,int>> edges, const std::map<std::pair<int,int>, int> &bond_orders);

// Functional Groups

void Parse_FG_to_Adj_Mat(std::vector<std::string> pdb_lines, AdjMat_int &tmp_adj_mat, std::vector<int> &atomic_numbers);
void Parse_PDBFile(std::string filename, std::vector<std::string> &atom_names, AdjMat_int &adjacency_matrix, std::vector<int> &atomic_numbers);
void GenerateKnownFGAdjMats(std::map<std::string, AdjMat_int> &known_fg_adj_mats, std::map<std::string, std::vector<int>> &known_fg_atomic_number_sets);
std::vector<std::string> IdentifyFunctionalGroups(std::string filename);
#endif
