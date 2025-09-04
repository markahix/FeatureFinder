#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include <vector>
#include <string>
// #include <iostream>
// #include <fstream>
// #include <map>
// #include <experimental/filesystem>
// #include <sstream>
#include <cmath>
#include <algorithm>
// #include <cctype>
// #include <locale>
// #include <cstdio>
// #include <memory>
// #include <stdexcept>
// #include <array>
// #include <iomanip>
// #include <ctime>
// #include <set>
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
    std::vector<int> atoms;
};

// periodic table information, mapped to elemental symbol?
std::vector<std::string> atomicSymbols;
int get_atomic_number(std::string symbol);
std::vector<double> vanDerWaalsRadii;


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
std::vector<Bond> GetBondsFromAtoms(std::vector<Atom> &atoms);
std::vector<Angle> GetAnglesFromBonds(std::vector<Atom> &atoms, std::vector<Bond> &bonds);
std::vector<Dihedral> GetDihedralsFromAngles(std::vector<Atom> &atoms, std::vector<Bond> &bonds, std::vector<Angle> &angles);
std::vector<Torsion> GetTorsionsFromAngles(std::vector<Atom> &atoms, std::vector<Bond> &bonds, std::vector<Angle> &angles);

// Pure Math Functions (using Points or vectors)
std::vector<double> cross_product(std::vector<double> vec_ab, std::vector<double> vec_bc);
std::vector<double> unit_vector(std::vector<double> vec_ab);
double dot_product(std::vector<double>vec_ab, std::vector<double>vec_bc);
double calc_distance(std::vector<double> vec_ab);
double calc_distance(Point A, Point B);
double calc_angle(Point A, Point B, Point C);
double calc_torsion(Point A, Point B, Point C, Point D);
#endif