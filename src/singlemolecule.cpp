/*
 * This portion handles all single-molecule processing such as identification of known functional groups and scaffolds.
 */

#include "singlemolecule.h"
#include "knownfgs.h"
#include "chemistry.h"

// Utilities
void trim(std::string &s)
{
    // rtrim(s);
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
    // ltrim(s);
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

bool val_in_array(int val, std::vector<int> my_array)
{
    for (int a : my_array)
    {
        if (a == val)
        {
            return true;
        }
    }
    return false;
}

bool Does_Vector_Contain_Vector(std::vector<int> bigger, std::vector<int> smaller)
{
    std::set<int> smaller_set(smaller.begin(), smaller.end());
    for (int s : smaller_set)
    {
        if (std::count(bigger.begin(), bigger.end(), s) < std::count(smaller.begin(), smaller.end(), s))
        {
            return false;
        }
    }
    return true;
}

bool any_val_in_array(std::set<int> vals, std::set<int>my_array)
{
    for (int v: vals)
    {
        for (int m: my_array)
        {
            if (m == v)
            {
                return true;
            }
        }
    }
    return false;
}

// Geometry
xyzCoords vector_between_points(xyzCoords pointA, xyzCoords pointB)
{
    xyzCoords vector;
    vector.x = pointB.x - pointA.x;
    vector.y = pointB.y - pointA.y;
    vector.z = pointB.z - pointA.z;
    return vector;
}

double norm_of_vector(xyzCoords vector)
{
    return sqrt(pow(vector.x,2) + pow(vector.y,2) + pow(vector.z,2));
}

double calc_distance_between_points(xyzCoords pointA, xyzCoords pointB)
{
    return norm_of_vector(vector_between_points(pointA, pointB));
}

// Adjacency Matrices
AdjacencyMatrix::AdjacencyMatrix(std::vector<int> incoming_nodes, std::vector<std::pair<int,int>> incoming_edges)
{
    node_weights = incoming_nodes;
    edges = incoming_edges;
    adjacency_matrix = {};
    std::vector<int> n_connections={};
    for (unsigned int i = 0; i < incoming_nodes.size(); i++)
    {
        n_connections.push_back(0);
    }

    for (std::pair<int,int> ed : edges)
    {
        n_connections[ed.first]++;
        n_connections[ed.second]++;
    }

    // regenerate atomic numbers by multiplying each value by 10 and then adding the n_connections value
    for (unsigned int i = 0; i < incoming_nodes.size(); i++)
    {
        incoming_nodes[i] = incoming_nodes[i]*10 + n_connections[i];
    }

    for (unsigned int i = 0; i < incoming_nodes.size(); i++)
    {
        adjacency_matrix.push_back({});
        for (unsigned int j = 0; j < incoming_nodes.size(); j++)
        {
            adjacency_matrix.back().push_back(0);
        }
    }

    // populate adjacency matrix
    for (std::pair<int,int> ed : edges)
    {
        int at1 = incoming_nodes[ed.first];
        int at2 = incoming_nodes[ed.second];
        int full_val = 0;
        if (at1 < at2)
        {
            full_val = 10000*at2 + at1;
        }
        else
        {
            full_val = 10000*at1 + at2;
        }
        adjacency_matrix[ed.first][ed.second] = full_val;
        adjacency_matrix[ed.second][ed.first] = full_val;
    }
}
AdjacencyMatrix::~AdjacencyMatrix()
{
    for (unsigned int i = 0; i < adjacency_matrix.size(); i++)
    {
        adjacency_matrix[i].clear();
    }
    adjacency_matrix.clear();
    node_weights.clear();
    edges.clear();
}
std::string AdjacencyMatrix::PrintAdjacencyMatrix()
{
    int cell_width = 8;
    std::stringstream buffer;
    buffer.str("");
    for (unsigned int i = 0; i < adjacency_matrix.size(); i++)
    {
        for (unsigned int j = 0; j < adjacency_matrix.size(); j++)
        {
            buffer << std::right << std::setw(cell_width) << std::setfill(' ') << adjacency_matrix[i][j] << " ";
        }
        buffer << std::endl;
    }
    return buffer.str();
}

namespace
{
    std::pair<int,int> normalized_edge(int a, int b)
    {
        if (a < b)
        {
            return {a, b};
        }
        return {b, a};
    }

    int infer_bond_order(int z1, int z2, double dist)
    {
        int a = std::min(z1, z2);
        int b = std::max(z1, z2);

        // These are intentionally simple, dependency-free heuristics.
        // They are designed to separate common single / double / aromatic-like
        // bonds well enough for functional-group matching without changing the
        // overall I/O scheme of the original code.
        if (a == 6 && b == 6)
        {
            if (dist <= 1.445) return 4; // aromatic-like C-C (captures benzene better than cyclohexane)
            return 1;
        }
        if (a == 6 && b == 7)
        {
            if (dist <= 1.315) return 2;
            if (dist <= 1.405) return 4;
            return 1;
        }
        if (a == 6 && b == 8)
        {
            if (dist <= 1.315) return 2;
            return 1;
        }
        if (a == 7 && b == 8)
        {
            if (dist <= 1.25) return 2;
            return 1;
        }
        if (a == 8 && b == 15)
        {
            if (dist <= 1.53) return 2;
            return 1;
        }
        if (a == 8 && b == 16)
        {
            if (dist <= 1.53) return 2;
            return 1;
        }
        return 1;
    }
}


/*
 * UPDATED IMPLEMENTATION NOTE
 * ---------------------------
 * The old Generate_Adjacency_Matrix() removed some atoms (notably legacy R-group
 * placeholders) while still indexing the matrix using the original atom indices.
 * That made the matrix/index relationship fragile and also made it hard to do a
 * correct graph search later.
 *
 * The updated version keeps a 1:1 relationship between:
 *   - atomic_numbers[i]
 *   - row/column i in the adjacency matrix
 *
 * We still preserve the legacy encoded node-label scheme:
 *   encoded_label = atomic_number * 10 + degree
 *
 * This keeps the existing I/O contract intact, while making the data safe for the
 * new matcher in main.cpp.
 *
 * Bond-type support can be added here later by changing the edge encoding. Right
 * now we still encode simple connectivity only.
 */
AdjMat_int Generate_Adjacency_Matrix(std::vector<int> &atomic_numbers, std::vector<std::pair<int,int>> edges)
{
    std::map<std::pair<int,int>, int> default_bond_orders = {};
    return Generate_Adjacency_Matrix(atomic_numbers, edges, default_bond_orders);
}

AdjMat_int Generate_Adjacency_Matrix(std::vector<int> &atomic_numbers, std::vector<std::pair<int,int>> edges, const std::map<std::pair<int,int>, int> &bond_orders)
{
    std::vector<int> n_connections(atomic_numbers.size(), 0);
    for (std::pair<int,int> ed : edges)
    {
        n_connections[ed.first]++;
        n_connections[ed.second]++;
    }

    for (unsigned int i = 0; i < atomic_numbers.size(); i++)
    {
        atomic_numbers[i] = atomic_numbers[i] * 10 + n_connections[i];
    }

    AdjMat_int adj_mat(atomic_numbers.size(), std::vector<int>(atomic_numbers.size(), 0));

    for (std::pair<int,int> ed : edges)
    {
        auto ek = normalized_edge(ed.first, ed.second);
        int bond_order = 1;
        auto it = bond_orders.find(ek);
        if (it != bond_orders.end())
        {
            bond_order = it->second;
        }

        // Adjacency entries now store bond-order labels directly.
        // Node labels still live in atomic_numbers[i] = atomic_number*10 + degree.
        adj_mat[ed.first][ed.second] = bond_order;
        adj_mat[ed.second][ed.first] = bond_order;
    }

    return adj_mat;
}



/*
 * UPDATED IMPLEMENTATION NOTE
 * ---------------------------
 * The parsing functions keep the same signatures, but the surrounding comments
 * now make the behavior explicit:
 *   - parse atoms from pseudo-PDB input
 *   - infer bonds from covalent-radius distance thresholds
 *   - build the legacy adjacency matrix + encoded node labels
 *
 * We intentionally leave the I/O contract untouched so the caller code can stay
 * almost identical.
 */
void Parse_FG_to_Adj_Mat(std::vector<std::string> pdb_lines, AdjMat_int &tmp_adj_mat, std::vector<int> &atomic_numbers)
{
    std::vector<xyzCoords> atomic_coords = {};
    std::vector<std::pair<int,int>> edges = {};
    std::map<std::pair<int,int>, int> bond_orders = {};

    for (unsigned int i = 0; i < pdb_lines.size(); i++)
    {
        std::string line = pdb_lines[i];
        atomic_coords.push_back(xyzCoords{stod(line.substr(30,8)), stod(line.substr(38,8)), stod(line.substr(46,8))});

        int element = get_atomic_number(line.substr(77,2));
        atomic_numbers.push_back(element);

        for (unsigned int j = 0; j < i; j++)
        {
            double curr_dist = calc_distance_between_points(atomic_coords[j], atomic_coords.back());
            double max_cov_dist = get_covalent_radius(atomic_numbers[j]) + get_covalent_radius(atomic_numbers.back());
            if (curr_dist < max_cov_dist)
            {
                edges.push_back({static_cast<int>(j), static_cast<int>(i)});
                bond_orders[normalized_edge(static_cast<int>(j), static_cast<int>(i))] = infer_bond_order(atomic_numbers[j], atomic_numbers.back(), curr_dist);
            }
        }
    }

    tmp_adj_mat = Generate_Adjacency_Matrix(atomic_numbers, edges, bond_orders);
}

void Parse_PDBFile(std::string filename, std::vector<std::string> &atom_names, AdjMat_int &adjacency_matrix, std::vector<int> &atomic_numbers)
{
    std::vector<std::string> pdb_lines = {};
    std::string line;
    std::ifstream infile(filename, std::ios::in);
    while (getline(infile,line)) // Parse PDB into atoms
    {
        if((line.substr(0,6)!= "ATOM  ") && (line.substr(0,6)!= "HETATM"))
        {
            continue;
        }
        pdb_lines.push_back(line);
    }
    std::vector<xyzCoords> atomic_coords = {};
    std::vector<std::pair<int,int>> edges = {};
    std::map<std::pair<int,int>, int> bond_orders = {};

    for (unsigned int i = 0; i < pdb_lines.size(); i++)
    {
        std::string line = pdb_lines[i];
        atomic_coords.push_back(xyzCoords{stod(line.substr(30,8)), stod(line.substr(38,8)), stod(line.substr(46,8))});

        std::string aname = line.substr(12,4);
        trim(aname);
        atom_names.push_back(aname);
        atomic_numbers.push_back(get_atomic_number(line.substr(77,2)));

        for (unsigned int j = 0; j < i; j++)
        {
            double curr_dist = calc_distance_between_points(atomic_coords[j], atomic_coords.back());
            double max_cov_dist = get_covalent_radius(atomic_numbers[j]) + get_covalent_radius(atomic_numbers.back());
            if (curr_dist < max_cov_dist)
            {
                edges.push_back({static_cast<int>(j), static_cast<int>(i)});
                bond_orders[normalized_edge(static_cast<int>(j), static_cast<int>(i))] = infer_bond_order(atomic_numbers[j], atomic_numbers.back(), curr_dist);
            }
        }
    }

    adjacency_matrix = Generate_Adjacency_Matrix(atomic_numbers, edges, bond_orders);
}


/*
 * UPDATED IMPLEMENTATION NOTE
 * ---------------------------
 * The linker error came from this definition being missing in the updated file.
 * The declaration remained in baseclasses.h and main.cpp still calls it, so the
 * linker quite correctly complained about an undefined reference.
 *
 * This implementation keeps the original function signature and behavior:
 *   - read the built-in pseudo-PDB functional-group blocks
 *   - convert each block to the legacy encoded adjacency-matrix form
 *   - store both the matrix and its encoded atom-label list in the provided maps
 *
 * We preserve the original size-prefixed key trick so larger functional groups
 * sort earlier and therefore override smaller ones later in the greedy pass.
 */
void GenerateKnownFGAdjMats(std::map<std::string, AdjMat_int> &known_fg_adj_mats,
                            std::map<std::string, std::vector<int>> &known_fg_atomic_number_sets)
{
    std::string fg_name;
    std::vector<std::string> block = {};

    for (const std::string &line : known_fgs)
    {
        if (line.find("REMARK") != std::string::npos)
        {
            fg_name = line.substr(7, line.size() - 7);
            trim(fg_name);
            continue;
        }

        if (line.empty())
        {
            if (block.empty())
            {
                continue;
            }

            std::vector<int> atomic_numbers = {};
            AdjMat_int tmp_adj_mat = {};
            Parse_FG_to_Adj_Mat(block, tmp_adj_mat, atomic_numbers);

            std::stringstream tmp_name;
            tmp_name.str("");
            tmp_name << std::setw(4) << std::setfill('0') << 100 - static_cast<int>(tmp_adj_mat.size()) << fg_name;

            known_fg_atomic_number_sets[tmp_name.str()] = atomic_numbers;
            known_fg_adj_mats[tmp_name.str()] = tmp_adj_mat;
            block.clear();
            continue;
        }

        if (line.find("ATOM") != std::string::npos)
        {
            block.push_back(line);
        }
    }
}

/*
 * UPDATED IMPLEMENTATION NOTE
 * ---------------------------
 * Legacy R-group placeholders are encoded as values < 10 because their atomic
 * number is 0 and only the degree contributes to the encoded label.
 *
 * The old code only skipped the exact value 1 in some places, which missed
 * non-terminal R groups. The updated screening step skips all encoded R-group
 * labels (< 10), which is consistent with the rest of the code.
 */

// Chemistry
// const std::vector<std::string> atomicSymbols = { "RR", // keeps symbol index aligned with atomic numbers.
//     "H" ,                                                                                                 "HE", 
//     "LI", "BE",                                                             "B" , "C" , "N" , "O" , "F" , "NE",
//     "NA", "MG",                                                             "AL", "SI", "P" , "S" , "CL", "AR", 
//     "K" , "CA", "SC", "TI", "V" , "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR", 
//     "RB", "SR", "Y" , "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I" , "XE", 
//     "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W" , "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN", 
//     "FR", "RA", "AC", "TH", "PA", "U" , "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN", "NH", "FL", "MC", "LV", "TS", "OG"
// };
// const std::vector<double> atomic_masses = {0.000, // NULL SPACE to align atomic number with mass index.
//       1.008,   4.003,   6.940,   9.012,  10.810,  12.011,  14.007,  15.999,  18.998,  20.180, // 1-10
//      22.990,  24.305,  26.982,  28.085,  30.974,  32.060,  35.450,  39.948,  39.098,  40.078, // 11-20
//      44.956,  47.867,  50.942,  51.996,  54.938,  55.845,  58.933,  58.693,  63.546,  65.380, // 21-30
//      69.723,  72.630,  74.922,  78.971,  79.904,  83.798,  85.468,  87.620,  88.906,  91.224, // 31-40
//      92.906,  95.950,  98.000, 101.070, 102.906, 106.420, 107.868, 112.414, 114.818, 118.710, // 41-50
//     121.760, 127.600, 126.904, 131.293, 132.905, 137.327, 138.905, 140.116, 140.908, 144.242, // 51-60
//     145.000, 150.360, 151.964, 157.250, 158.925, 162.500, 164.930, 167.259, 168.934, 173.045, // 61-70
//     174.967, 178.490, 180.948, 183.840, 186.207, 190.230, 192.217, 195.084, 196.967, 200.592, // 71-80
//     204.380, 207.200, 208.980, 209.000, 210.000, 222.000, 223.000, 226.000, 227.000, 232.038, // 81-90
//     231.036, 238.029, 237.000, 244.000, 243.000, 247.000, 247.000, 251.000, 252.000, 257.000, // 91-100
//     258.000, 259.000, 266.000, 267.000, 268.000, 269.000, 270.000, 269.000, 278.000, 281.000, // 100-110
//     282.000, 285.000, 286.000, 289.000, 290.000, 293.000, 294.000, 294.000                    // 111-118
// };
// const std::vector<double> covalent_radii = {  0.76, // included 0.00 element such that array indices align with atomic numbers.
//     0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,
//     1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76,
//     1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.21, 1.38,
//     1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75,
//     1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.41, 1.38,
//     1.35, 1.33, 1.33, 1.31, 2.44, 2.23, 2.01, 1.99, 1.98, 1.98,
//     1.96, 1.95, 1.93, 1.92, 1.92, 1.89, 1.90, 1.87, 1.86, 1.85,
//     1.85, 1.80, 1.78, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36,
//     1.33, 1.31, 1.29, 1.28, 1.27, 1.25, 2.60, 2.53, 2.21, 2.15,
//     2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.80, 1.77, 1.74, 1.73,
//     1.72, 1.68, 1.65, 1.57, 1.49, 1.43, 1.41, 1.34, 1.29, 1.28,
//     1.21, 1.22, 1.20, 1.20, 1.20, 1.20, 1.20
// };
// const std::vector<double> vanDerWaalsRadii = {0.00,
//     1.20, 1.40, 1.82, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54, // 1–10
//     2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88, 2.75, 2.31, // 11–20
//     2.11, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 1.63, 1.40, 1.39, // 21–30
//     1.87, 2.11, 1.85, 1.90, 1.85, 2.02, 3.03, 2.49, 2.00, 2.06, // 31–40
//     1.98, 2.17, 2.16, 2.06, 2.02, 2.05, 2.03, 2.20, 2.20, 2.25, // 41–50
//     2.20, 2.06, 1.98, 2.16, 3.43, 2.68, 2.40, 2.35, 2.39, 2.29, // 51–60
//     2.36, 2.29, 2.33, 2.37, 2.21, 2.29, 2.16, 2.35, 2.27, 2.32, // 61–70
//     2.31, 2.25, 2.38, 2.32, 2.33, 2.30, 2.26, 2.30, 2.32, 2.26, // 71–80
//     1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 3.48, 2.83, 2.00, 2.40, // 81–90
//     2.00, 1.86, 1.90, 1.87, 1.80, 1.69, 1.68, 1.68, 1.65, 1.67, // 91–100
//     1.73, 1.76, 1.61, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  // 101–110
//     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0                // 111–118
// };
// int get_atomic_number(std::string symbol)
// {
//     return std::distance(atomicSymbols.begin(), std::find(atomicSymbols.begin(), atomicSymbols.end(), symbol));
// }
// double get_covalent_radius(std::string symbol)
// {
//     return covalent_radii[get_atomic_number(symbol)]*1.05;
// }
// double get_covalent_radius(int number)
// {
//     return covalent_radii[number]*1.05;
// }
// double get_vdw_radius(std::string symbol)
// {
//     return vanDerWaalsRadii[get_atomic_number(symbol)]*1.05;
// }
// double get_vdw_radius(int number)
// {
//     return vanDerWaalsRadii[number]*1.05;
// }
// double get_atomic_mass(std::string symbol)
// {
//     return atomic_masses[get_atomic_number(symbol)];
// }
// double get_atomic_mass(int number)
// {
//     return atomic_masses[number];
// }


inline bool is_r_group_label(int encoded_label)
{
    // Legacy convention: RR has atomic number 0, so after encoding the value is just its degree (0..9).
    return encoded_label < 10;
}

inline int encoded_atomic_number(int encoded_label)
{
    return encoded_label / 10;
}

SearchGraph build_search_graph(const AdjMat_int &adj_mat, const std::vector<int> &encoded_labels)
{
    SearchGraph g;
    g.labels = encoded_labels;
    g.neighbors.assign(adj_mat.size(), {});
    g.edge_label.assign(adj_mat.size(), {});
    g.degree.assign(adj_mat.size(), 0);

    for (unsigned int i = 0; i < adj_mat.size(); i++)
    {
        for (unsigned int j = i + 1; j < adj_mat[i].size(); j++)
        {
            if (adj_mat[i][j] != 0)
            {
                g.neighbors[i].push_back(static_cast<int>(j));
                g.neighbors[j].push_back(static_cast<int>(i));
                g.edge_label[i][static_cast<int>(j)] = adj_mat[i][j];
                g.edge_label[j][static_cast<int>(i)] = adj_mat[i][j];
                g.degree[i]++;
                g.degree[j]++;
            }
        }
    }
    return g;
}

PreparedQuery prepare_query_graph(const AdjMat_int &fgam, const std::vector<int> &fg_atom_list)
{
    PreparedQuery pq;

    std::vector<int> explicit_orig_indices = {};
    std::vector<int> orig_to_explicit(fg_atom_list.size(), -1);

    for (unsigned int q = 0; q < fg_atom_list.size(); q++)
    {
        if (!is_r_group_label(fg_atom_list[q]))
        {
            orig_to_explicit[q] = static_cast<int>(explicit_orig_indices.size());
            explicit_orig_indices.push_back(static_cast<int>(q));
            pq.encoded_label.push_back(fg_atom_list[q]);
        }
    }

    pq.orig_index_of_explicit = explicit_orig_indices;
    pq.neighbors.assign(explicit_orig_indices.size(), {});
    pq.edge_label.assign(explicit_orig_indices.size(), {});
    pq.required_external_neighbors.assign(explicit_orig_indices.size(), 0);
    pq.total_required_degree.assign(explicit_orig_indices.size(), 0);

    for (unsigned int qe = 0; qe < explicit_orig_indices.size(); qe++)
    {
        int q_orig = explicit_orig_indices[qe];
        for (unsigned int q2 = 0; q2 < fgam[q_orig].size(); q2++)
        {
            if (fgam[q_orig][q2] == 0)
            {
                continue;
            }

            if (is_r_group_label(fg_atom_list[q2]))
            {
                // R-group placeholders act like wildcard external neighbors.
                pq.required_external_neighbors[qe]++;
            }
            else
            {
                int qe2 = orig_to_explicit[q2];
                if (qe2 >= 0)
                {
                    pq.neighbors[qe].push_back(qe2);
                    pq.edge_label[qe][qe2] = fgam[q_orig][q2];
                }
            }
        }

        pq.total_required_degree[qe] = static_cast<int>(pq.neighbors[qe].size()) + pq.required_external_neighbors[qe];
    }

    for (unsigned int qe = 0; qe < pq.neighbors.size(); qe++)
    {
        std::sort(pq.neighbors[qe].begin(), pq.neighbors[qe].end());
        pq.neighbors[qe].erase(std::unique(pq.neighbors[qe].begin(), pq.neighbors[qe].end()), pq.neighbors[qe].end());
    }

    return pq;
}

std::map<int,int> neighbor_atomic_histogram(const SearchGraph &g, int v)
{
    std::map<int,int> hist = {};
    for (int n : g.neighbors[v])
    {
        hist[encoded_atomic_number(g.labels[n])]++;
    }
    return hist;
}

std::map<int,int> query_neighbor_atomic_histogram(const PreparedQuery &qg, int qv)
{
    std::map<int,int> hist = {};
    for (int qn : qg.neighbors[qv])
    {
        hist[encoded_atomic_number(qg.encoded_label[qn])]++;
    }
    return hist;
}

bool histogram_subset(const std::map<int,int> &need, const std::map<int,int> &have)
{
    for (auto const &kv : need)
    {
        auto it = have.find(kv.first);
        int cnt = (it == have.end()) ? 0 : it->second;
        if (cnt < kv.second)
        {
            return false;
        }
    }
    return true;
}

bool target_node_compatible(int q_label, int q_required_degree, const std::map<int,int> &q_hist,
                            const SearchGraph &target, int tv, const std::set<int> &declined)
{
    if (declined.find(tv) != declined.end())
    {
        return false;
    }

    if (is_r_group_label(q_label))
    {
        return false; // explicit matcher never maps query R-group placeholders directly.
    }

    if (encoded_atomic_number(q_label) != encoded_atomic_number(target.labels[tv]))
    {
        return false;
    }

    // Require exact degree equality.
    // This prevents six sp3 carbons in cyclohexane from matching six sp2 carbons
    // in benzene, and more generally respects the full local valence implied by
    // the query plus any wildcard R-group attachments.
    if (target.degree[tv] != q_required_degree)
    {
        return false;
    }

    std::map<int,int> t_hist = neighbor_atomic_histogram(target, tv);
    if (!histogram_subset(q_hist, t_hist))
    {
        return false;
    }

    return true;
}

std::vector<int> location_of_fg_in_mol(AdjMat_int fgam, AdjMat_int molam, std::vector<int> fg_atom_list, std::vector<int> mol_atom_list, std::vector<int> already_identified_indices)
{
    /*
     * UPDATED IMPLEMENTATION NOTE
     * ---------------------------
     * Same signature, new algorithm~
     *
     * Old behavior:
     *   - enumerate every k-subset of remaining molecule atoms
     *   - build a submatrix for each subset
     *   - compare the whole submatrix to the FG matrix
     *
     * New behavior:
     *   - build sparse graph views for query and target
     *   - ignore R-group placeholders as explicit match targets
     *   - convert each R-group into a required “external neighbor count” constraint
     *   - run recursive subgraph matching with pruning
     *
     * The output contract stays the same: return the molecule indices for one found FG
     * occurrence, or {} if nothing matches.
     */

    SearchGraph target = build_search_graph(molam, mol_atom_list);
    PreparedQuery query = prepare_query_graph(fgam, fg_atom_list);

    if (query.orig_index_of_explicit.empty())
    {
        return {};
    }

    std::set<int> declined(already_identified_indices.begin(), already_identified_indices.end());

    // Precompute lightweight local signatures for the explicit query atoms.
    std::vector<std::map<int,int>> query_neighbor_hists = {};
    for (unsigned int q = 0; q < query.orig_index_of_explicit.size(); q++)
    {
        query_neighbor_hists.push_back(query_neighbor_atomic_histogram(query, static_cast<int>(q)));
    }

    // Initial candidate domains: only plausible target atoms are kept.
    std::vector<std::vector<int>> domains(query.orig_index_of_explicit.size(), std::vector<int>{});
    for (unsigned int q = 0; q < query.orig_index_of_explicit.size(); q++)
    {
        for (int tv = 0; tv < static_cast<int>(target.labels.size()); tv++)
        {
            if (target_node_compatible(query.encoded_label[q], query.total_required_degree[q], query_neighbor_hists[q], target, tv, declined))
            {
                domains[q].push_back(tv);
            }
        }

        if (domains[q].empty())
        {
            return {};
        }
    }

    // Query-order heuristic: start with the most constrained atom first.
    std::vector<int> order(query.orig_index_of_explicit.size(), 0);
    for (unsigned int i = 0; i < order.size(); i++)
    {
        order[i] = static_cast<int>(i);
    }
    std::sort(order.begin(), order.end(), [&](int a, int b)
    {
        if (domains[a].size() != domains[b].size())
        {
            return domains[a].size() < domains[b].size();
        }
        return query.total_required_degree[a] > query.total_required_degree[b];
    });

    std::vector<int> mapping(query.orig_index_of_explicit.size(), -1); // explicit query vertex -> target vertex
    std::set<int> used_targets = {};
    std::vector<int> found_result = {};

    std::function<bool(int)> backtrack = [&](int depth) -> bool
    {
        if (depth == static_cast<int>(order.size()))
        {
            // Final check for wildcard R-groups: each explicit query atom that had one or more
            // R-group placeholders must have enough neighbors outside the mapped explicit query core.
            for (unsigned int q = 0; q < query.orig_index_of_explicit.size(); q++)
            {
                int tv = mapping[q];
                int mapped_explicit_neighbors = 0;
                for (int qn : query.neighbors[q])
                {
                    int tn = mapping[qn];
                    for (int real_n : target.neighbors[tv])
                    {
                        if (real_n == tn)
                        {
                            mapped_explicit_neighbors++;
                            break;
                        }
                    }
                }

                int available_external = target.degree[tv] - mapped_explicit_neighbors;
                if (available_external < query.required_external_neighbors[q])
                {
                    return false;
                }
            }

            found_result.clear();
            for (int tv : mapping)
            {
                found_result.push_back(tv);
            }
            std::sort(found_result.begin(), found_result.end());
            found_result.erase(std::unique(found_result.begin(), found_result.end()), found_result.end());
            return true;
        }

        int qv = order[depth];

        for (int tv : domains[qv])
        {
            if (used_targets.find(tv) != used_targets.end())
            {
                continue;
            }

            bool consistent = true;

            // Incremental edge consistency: whenever a query neighbor is already mapped,
            // the corresponding target edge must exist.
            for (int qn : query.neighbors[qv])
            {
                if (mapping[qn] == -1)
                {
                    continue;
                }

                int tn = mapping[qn];
                auto it = target.edge_label[tv].find(tn);
                if (it == target.edge_label[tv].end())
                {
                    consistent = false;
                    break;
                }

                int required_bond_label = query.edge_label[qv].at(qn);
                if (it->second != required_bond_label)
                {
                    consistent = false;
                    break;
                }
            }

            if (!consistent)
            {
                continue;
            }

            mapping[qv] = tv;
            used_targets.insert(tv);

            if (backtrack(depth + 1))
            {
                return true;
            }

            used_targets.erase(tv);
            mapping[qv] = -1;
        }

        return false;
    };

    if (backtrack(0))
    {
        return found_result;
    }

    return {};
}
std::string get_residue_name(std::string filename)
{
    std::string line;
    std::string resname="LIG";
    std::ifstream ifile(filename,std::ios::in);
    while (getline(ifile,line))
    {
        if ((line.find("ATOM") != std::string::npos) || (line.find("HETATM") != std::string::npos))
        {
            resname = line.substr(17,3);
            break;
        }
    }
    ifile.close();
    return resname;
}

std::vector<std::string> IdentifyFunctionalGroups(std::string filename)
{   
    std::vector<std::string> ReturningAGIMUSLines = {};
    std::string residue_name = get_residue_name(filename);
    // Initialize Variables
    std::map<std::string, AdjMat_int> known_fg_adj_mats = {};
    std::map<std::string, std::vector<int>> known_fg_atomic_number_sets = {};
    std::vector<std::string> mol_atom_names = {};
    std::vector<int> mol_atomic_numbers = {};
    AdjMat_int mol_adj_mat;

    GenerateKnownFGAdjMats(known_fg_adj_mats, known_fg_atomic_number_sets);
    std::cout << "Currently " << known_fg_adj_mats.size() << " known functional groups." << std::endl;

    Parse_PDBFile(filename, mol_atom_names, mol_adj_mat, mol_atomic_numbers);
    std::cout << "Loaded " << filename << std::endl;
    std::cout << "Identified " << mol_atom_names.size() << " atoms." << std::endl;

    // std::vector<std::string> possible_fgs_to_check = get_possible_fgs(known_fg_atomic_number_sets, mol_atomic_numbers);


    // Keep the same overall structure as the original code: prepare per-FG atom options,
    // then repeatedly search and remove identified atoms from consideration.
    std::map<std::string, std::vector<int>> atoms_by_fg_option = {};
    for (std::map<std::string, AdjMat_int>::iterator fg_iter=known_fg_adj_mats.begin(); fg_iter != known_fg_adj_mats.end(); fg_iter++)
    // for (std::string poss_fg : possible_fgs_to_check)
    {
        std::string poss_fg = fg_iter->first;
        atoms_by_fg_option[poss_fg] = {};
        for (int i = 0; i < static_cast<int>(mol_atomic_numbers.size()); i++)
        {
            for (int qlab : known_fg_atomic_number_sets[poss_fg])
            {
                if (is_r_group_label(qlab))
                {
                    continue;
                }
                if (mol_atomic_numbers[i] == qlab)
                {
                    atoms_by_fg_option[poss_fg].push_back(i);
                    break;
                }
            }
        }
    }
    std::cout << "Filled " << atoms_by_fg_option.size() << " atoms_by_fg_option entries"<<std::endl;

    std::vector<int> already_identified_indices = {};
    std::set<int> claimed_atoms = {};

    /*
     * SUMMARY STORAGE
     * ---------------
     * We keep the original streaming output above, but also collect a compact
     * per-functional-group summary so the final result is much easier to read.
     *
     * summary_atom_names[fg_name] = list of matches, each match = list of atom names
     * summary_atom_indices[fg_name] = list of matches, each match = list of atom indices
     *
     * This is intentionally additive only: it does NOT change the search logic,
     * matching behavior, or I/O contract. It only records what was already found.
     */
    std::map<std::string, std::vector<std::vector<std::string>>> summary_atom_names = {};
    std::map<std::string, std::vector<std::vector<int>>> summary_atom_indices = {};

    for (std::map<std::string, AdjMat_int>::iterator fg_iter=known_fg_adj_mats.begin(); fg_iter != known_fg_adj_mats.end(); fg_iter++)
    // for (std::string poss_fg : possible_fgs_to_check)
    {
        std::string poss_fg = fg_iter->first;
        std::vector<int> tmp_declined_indices = {};

        for (int i = 0; i < static_cast<int>(mol_atomic_numbers.size()); i++)
        {
            // Exclusion rule. Once an atom has been claimed by an earlier, larger
            // functional group, all later searches must skip it completely.
            if (claimed_atoms.find(i) != claimed_atoms.end())
            {
                tmp_declined_indices.push_back(i);
                continue;
            }
            if (!val_in_array(i, atoms_by_fg_option[poss_fg]))
            {
                tmp_declined_indices.push_back(i);
                continue;
            }
        }

        AdjMat_int curr_fg_adjmat = known_fg_adj_mats[poss_fg];
        std::vector<int> curr_fg_atom_list = known_fg_atomic_number_sets[poss_fg];

        std::vector<int> new_fg_indices;
        do
        {
            new_fg_indices = location_of_fg_in_mol(curr_fg_adjmat, mol_adj_mat, curr_fg_atom_list, mol_atomic_numbers, tmp_declined_indices);
            if (new_fg_indices.size() == 0)
            {
                continue;
            }

            /*
             * RECORD THIS MATCH FOR THE END-OF-RUN SUMMARY.
             * We store both atom names and raw atom indices so the user can
             * quickly see either representation later.
             */
            std::vector<std::string> matched_atom_names = {};
            std::vector<int> matched_atom_indices = {};

            bool overlaps_claimed = false;
            for (int a : new_fg_indices)
            {
                if (claimed_atoms.find(a) != claimed_atoms.end())
                {
                    overlaps_claimed = true;
                    break;
                }
            }
            if (overlaps_claimed)
            {
                new_fg_indices.clear();
                continue;
            }

            for (int a : new_fg_indices)
            {
                matched_atom_names.push_back(mol_atom_names[a]);
                matched_atom_indices.push_back(a);
                already_identified_indices.push_back(a);
                claimed_atoms.insert(a);
                tmp_declined_indices.push_back(a);
            }

            summary_atom_names[poss_fg].push_back(matched_atom_names);
            summary_atom_indices[poss_fg].push_back(matched_atom_indices);
        } while (new_fg_indices.size() > 0);
    }

    /*
     * FINAL HUMAN-READABLE SUMMARY
     * ----------------------------
     * Print one compact section per functional group showing:
     *   - how many matches were found
     *   - atom names for each match
     *   - atom indices for each match
     *
     * This makes the output much easier to scan than the streaming log above.
     */
    std::cout << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "FUNCTIONAL GROUP SUMMARY" << std::endl;
    std::cout << "============================================================" << std::endl;

    bool printed_any_summary = false;
    for (std::map<std::string, AdjMat_int>::iterator fg_iter=known_fg_adj_mats.begin(); fg_iter != known_fg_adj_mats.end(); fg_iter++)
    {
        // To Do:  Build in a std::vector<std::string> AGIMUS_FunctionalGroups lines which contain the 
        // name of the functional group, followed by an Amber-formatted atom mask of the selection.
        // use std::stringstream buffer to build each line.
        std::stringstream buffer;
        buffer.str("");
        std::string poss_fg = fg_iter->first;
        auto it_names = summary_atom_names.find(poss_fg);
        auto it_indices = summary_atom_indices.find(poss_fg);
        if (it_names == summary_atom_names.end() || it_indices == summary_atom_indices.end())
        {
            continue;
        }
        if (it_names->second.empty())
        {
            continue;
        }

        printed_any_summary = true;
        std::cout << poss_fg << ": " << it_names->second.size() << " match(es)" << std::endl;
        
        for (std::size_t m = 0; m < it_names->second.size(); m++)
        {
            std::cout << "  [" << (m + 1) << "] names:   ";
            buffer <<"AGIMUS_FUNCGRP " << poss_fg << ":" << residue_name<<"@";
            for (const std::string &nm : it_names->second[m])
            {
                std::cout << nm << " ";
                buffer << nm <<",";
            }
            std::cout << std::endl;
            ReturningAGIMUSLines.push_back(buffer.str());
            buffer.str("");

            std::cout << "      indices: ";
            for (int idx : it_indices->second[m])
            {
                std::cout << idx << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    if (!printed_any_summary)
    {
        std::cout << "No functional groups were found." << std::endl;
    }
    return ReturningAGIMUSLines;

}