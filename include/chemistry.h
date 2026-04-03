#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>

// #include <iostream>
// #include <experimental/filesystem>
// #include <cctype>
// #include <locale>
// #include <cstdio>
// #include <memory>
// #include <stdexcept>
// #include <array>

// #include <ctime>
#include <set>
// #include <cstdlib>


// Baseline Data Structures
struct Point
{
    double x;
    double y;
    double z;
};

struct Atom
{
    int index;
    double x;
    double y;
    double z;
    int atomic_number;
    std::map<std::string,std::string> info;
};

struct Bond
{
    Atom* atom_1;
    Atom* atom_2;
    double length;
    double bond_order;
};

struct Angle
{
    Atom* atom_1;
    Atom* atom_2;
    Atom* atom_3;
    double angle;
};

struct Torsion
{
    Atom* atom_1;
    Atom* atom_2;
    Atom* atom_3;
    Atom* atom_4;
    double angle;
};

struct Dihedral
{
    Atom* atom_1;
    Atom* atom_2;
    Atom* atom_3;
    Atom* atom_4;
    double angle;
};

struct Ring
{
    std::vector<Atom*> atoms;
};

class Molecule
{
    private:
        // Data Storage
        std::vector<Atom> atoms;
        std::vector<Bond> bonds;
        std::vector<Angle> angles;
        std::vector<Torsion> torsions;
        std::vector<Dihedral> dihedrals;
        std::vector<Ring> rings;

        // internal functions based on the addition of atoms.
        void _add_bond(Atom* A, Atom* B);
        void _add_angle(Atom* A, Atom* B, Atom* C);
        void _add_torsion(Atom* A, Atom* B, Atom* C, Atom* D);
        void _add_dihedral(Atom* A, Atom* B, Atom* C, Atom* D);
        void _add_ring(std::vector<Atom*> atoms);
        
        // value checks
        bool _is_aromatic(Ring R);
        bool _is_ring_atom(Atom* A);
        bool _is_halogen(Atom* A);
        
        // Filetype-specific parsers
        void _parse_from_pdb(std::string filename);
        void _parse_from_xyz(std::string filename);
        void _parse_from_mol2(std::string filename);
        void _parse_from_sdf(std::string filename);
        void _parse_from_prmtop_rst7(std::string filename1, std::string filename2);

    public:
        //Overloaded Constructors
        Molecule();
        Molecule(std::string filename);
        Molecule(std::string filename1,std::string filename2);
        Molecule(std::vector<std::string> mol_chunk);
        
        //Destructor
        ~Molecule();

        // Manually add atoms one-by-one or get called by constructor with arguments.
        void AddAtom(Atom &A);
        void SavePDB(std::string filename);
        void SaveXYZ(std::string filename);
        void SaveMOL2(std::string filename);
};

// periodic table information, mapped to elemental symbol?
int get_atomic_number(std::string symbol);

// molecular calculations
double calc_distance(double x1, double y1, double z1, double x2, double y2, double z2);
double calc_distance(std::vector<double> coords1, std::vector<double> coords2);
double calc_distance(Atom A, Atom B);
double calc_angle(Atom A, Atom B, Atom C);
double calc_torsion(Atom A, Atom B, Atom C, Atom D);
double calc_dihedral(Atom A, Atom B, Atom C, Atom D);
bool is_bonded(std::string element1,std::string element2, double x1, double y1, double z1, double x2, double y2, double z2);
bool is_bonded(Atom atom1, Atom atom2);

// Process Data Structures into larger Data Structures.
Point AtomToPoint(Atom A);
std::vector<Bond> GetBondsFromAtoms(std::vector<Atom> &atoms);
std::vector<Angle> GetAnglesFromBonds(std::vector<Atom> &atoms, std::vector<Bond> &bonds);
std::vector<Dihedral> GetDihedralsFromAngles(std::vector<Atom> &atoms, std::vector<Bond> &bonds, std::vector<Angle> &angles);
std::vector<Torsion> GetTorsionsFromAngles(std::vector<Atom> &atoms, std::vector<Bond> &bonds, std::vector<Angle> &angles);
std::vector<Ring> GetRingsFromAll(std::vector<Atom> &atoms, std::vector<Bond> &bonds, std::vector<Angle> &angles, std::vector<Torsion> &torsions);

// Pure Math Functions (using Points or vectors)
std::vector<double> cross_product(std::vector<double> vec_ab, std::vector<double> vec_bc);
std::vector<double> unit_vector(std::vector<double> vec_ab);
double dot_product(std::vector<double>vec_ab, std::vector<double>vec_bc);
double calc_distance(std::vector<double> vec_ab);
double calc_distance(Point A, Point B);
double calc_angle(Point A, Point B, Point C);
double calc_torsion(Point A, Point B, Point C, Point D);
#endif