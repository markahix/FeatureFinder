#include "chemistry.h"

const std::vector<std::string> atomicSymbols = { "NULL",
    "H",                                                                                                  "HE", 
    "LI", "BE",                                                             "B",  "C",  "N",  "O",  "F",  "NE",
    "NA", "MG",                                                             "AL", "SI", "P",  "S",  "CL", "AR", 
    "K",  "CA", "SC", "TI", "V",  "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR", 
    "RB", "SR",  "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I",  "XE", 
    "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W",  "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN", 
    "FR", "RA", "AC", "TH", "PA", "U",  "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN", "NH", "FL", "MC", "LV", "TS", "OG"
};

const std::vector<double> vanDerWaalsRadii = {0.00,
    1.20, 1.40, 1.82, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54, // 1–10
    2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88, 2.75, 2.31, // 11–20
    2.11, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 1.63, 1.40, 1.39, // 21–30
    1.87, 2.11, 1.85, 1.90, 1.85, 2.02, 3.03, 2.49, 2.00, 2.06, // 31–40
    1.98, 2.17, 2.16, 2.06, 2.02, 2.05, 2.03, 2.20, 2.20, 2.25, // 41–50
    2.20, 2.06, 1.98, 2.16, 3.43, 2.68, 2.40, 2.35, 2.39, 2.29, // 51–60
    2.36, 2.29, 2.33, 2.37, 2.21, 2.29, 2.16, 2.35, 2.27, 2.32, // 61–70
    2.31, 2.25, 2.38, 2.32, 2.33, 2.30, 2.26, 2.30, 2.32, 2.26, // 71–80
    1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 3.48, 2.83, 2.00, 2.40, // 81–90
    2.00, 1.86, 1.90, 1.87, 1.80, 1.69, 1.68, 1.68, 1.65, 1.67, // 91–100
    1.73, 1.76, 1.61, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  // 101–110
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0                // 111–118
};

#include <vector>

const std::vector<double> covalent_radii = {  0.00, // included 0.00 element such that array indices align with atomic numbers.
    0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,
    1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76,
    1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.21, 1.38,
    1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75,
    1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.41, 1.38,
    1.35, 1.33, 1.33, 1.31, 2.44, 2.23, 2.01, 1.99, 1.98, 1.98,
    1.96, 1.95, 1.93, 1.92, 1.92, 1.89, 1.90, 1.87, 1.86, 1.85,
    1.85, 1.80, 1.78, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36,
    1.33, 1.31, 1.29, 1.28, 1.27, 1.25, 2.60, 2.53, 2.21, 2.15,
    2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.80, 1.77, 1.74, 1.73,
    1.72, 1.68, 1.65, 1.57, 1.49, 1.43, 1.41, 1.34, 1.29, 1.28,
    1.21, 1.22, 1.20, 1.20, 1.20, 1.20, 1.20
};


int get_atomic_number(std::string symbol)
{
    for (unsigned int i = 0; i < atomicSymbols.size(); i++)
    {
        if (symbol == atomicSymbols[i])
        {
            return i;
        }
    }
    return 0;
}

double get_vdw_radius(int atomic_number)
{
    return vanDerWaalsRadii[atomic_number];
}
double get_vdw_radius(std::string symbol)
{
    return vanDerWaalsRadii[get_atomic_number(symbol)];
}
double get_covalent_radius(int atomic_number)
{
    return covalent_radii[atomic_number];
}
double get_covalent_radius(std::string symbol)
{
    return covalent_radii[get_atomic_number(symbol)];
}


double calc_distance(std::vector<double> vec_ab)
{
    double sum = 0;
    for (double val : vec_ab)
    {
        sum += pow(val,2);
    }
    return sqrt(sum);
}

double calc_distance(Point A, Point B)
{
    return calc_distance(A.x,A.y,A.z,B.x,B.y,B.z);
}

double calc_distance(Atom A, Atom B)
{
    return calc_distance(A.x,A.y,A.z,B.x,B.y,B.z);
}

double calc_distance(std::vector<double> coords1, std::vector<double> coords2)
{
    return calc_distance(coords1[0],coords1[1],coords1[2],coords2[0],coords2[1],coords2[2]);
}

double calc_distance( double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt( pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2) );
}

bool is_bonded(std::string element1, std::string element2, double x1, double y1, double z1, double x2, double y2, double z2)
{
    double vdw1 = get_vdw_radius(element1);
    double vdw2 = get_vdw_radius(element2);
    double max_dist = (vdw1 + vdw2)/2;
    double curr_dist = calc_distance(x1, y1, z1, x2, y2, z2);
    if (curr_dist <  max_dist)
    {
        return true;
    }
    return false;
}

bool is_bonded(Atom atom1, Atom atom2)
{
    double vdw1 = get_covalent_radius(atomicSymbols[atom1.atomic_number]);
    double vdw2 = get_covalent_radius(atomicSymbols[atom2.atomic_number]);
    double max_dist = (vdw1 + vdw2);
    double curr_dist = calc_distance(atom1,atom2);
    if (curr_dist <  max_dist)
    {
        return true;
    }
    return false;
}

std::vector<double> calc_vector(Point A, Point B)
{
    return {B.x-A.x, B.y-A.y, B.z-A.z};
}

double dot_product(std::vector<double>vec_ab, std::vector<double>vec_bc)
{
    double sum = 0;
    for (unsigned int i = 0; i < vec_ab.size(); i++)
    {
        sum += (vec_ab[i]*vec_bc[i]);
    }
    return sum;
}

double calc_angle(std::vector<double> vec_ab, std::vector<double> vec_bc)
{
    double len_ab = calc_distance(vec_ab);    // Magnitude (length) of vec_AB
    double len_bc = calc_distance(vec_bc);    // Magnitude (length) of vec_BC
    double ab_dot_bc =  dot_product(vec_ab,vec_bc);  // Dot product of vec_AB and vec_BC
    double theta_rad = acos((ab_dot_bc)/(len_ab * len_bc));    // theta = arccos(ab_dot_bc/(len_ab * len_bc))
    return theta_rad * 57.2958;    // convert to degrees and return (1 rad ~ 57.2958 degrees)
}

double calc_angle(Point A, Point B, Point C)
{
    std::vector<double> vec_ab = calc_vector(A, B);    // Vector from A to B
    std::vector<double> vec_bc = calc_vector(B, C);    // Vector from B to C
    return calc_angle(vec_ab,vec_bc);
}

Point AtomToPoint(Atom A)
{
    Point p={A.x,A.y,A.z};
    return p;
}

double calc_angle(Atom A, Atom B, Atom C)
{
    return calc_angle(AtomToPoint(A),AtomToPoint(B),AtomToPoint(C));
}

std::vector<double> cross_product(std::vector<double> vec_ab,std::vector<double> vec_bc)
{
    double i=(vec_ab[1]*vec_bc[2] - vec_ab[2]*vec_bc[1]);
    double j=(vec_ab[0]*vec_bc[2] - vec_ab[2]*vec_bc[0]);
    double k=(vec_ab[0]*vec_bc[1] - vec_ab[1]*vec_bc[0]);
    return {i,j,k};
}

std::vector<double> unit_vector(std::vector<double> vec_ab)
{
    double mag = calc_distance(vec_ab);
    std::vector<double> unit_vec = {};
    for (double val : vec_ab)
    {
        unit_vec.push_back(val/mag);
    }
    return unit_vec;
}

double calc_torsion(Point A, Point B, Point C, Point D)
{
    // obtain three vectors
    std::vector<double> vec_ab = calc_vector(A, B);    // Vector from A to B
    std::vector<double> vec_bc = calc_vector(B, C);    // Vector from B to C
    std::vector<double> vec_cd = calc_vector(C, D);    // Vector from C to D

    // obtain normal vectors for planes defined by vector pairs AB,BC and BC,CD
    std::vector<double> n1 = unit_vector(cross_product(vec_ab,vec_bc));
    std::vector<double> n2 = unit_vector(cross_product(vec_bc,vec_cd));
    std::vector<double> m1 = cross_product(n1,unit_vector(vec_bc));

    // obtain angle between normals
    double x = dot_product(n1,n2);
    double y = dot_product(m1,n2);
    return atan2(x, y);
}

double calc_torsion(Atom A, Atom B, Atom C, Atom D)
{
    Point pA={A.x,A.y,A.z};
    Point pB={B.x,B.y,B.z};
    Point pC={C.x,C.y,C.z};
    Point pD={D.x,D.y,D.z};
    return calc_torsion(pA,pB,pC,pD);
}

double calc_dihedral(Atom A, Atom B, Atom C, Atom D)
{
    Point pA={A.x,A.y,A.z};
    Point pB={B.x,B.y,B.z};
    Point pC={C.x,C.y,C.z};
    Point pD={D.x,D.y,D.z};
    return calc_torsion(pA,pB,pC,pD);
}

std::vector<Bond> GetBondsFromAtoms(std::vector<Atom> &atoms)
{
    std::vector<Bond> bonds = {};
    for (unsigned int i = 0; i < atoms.size()-1; i++ )
    {
        for (unsigned int j = i+1; j < atoms.size(); j++ )
        {
            Atom* atom1 = &atoms[i];
            Atom* atom2 = &atoms[j];
            std::string element1 = atomicSymbols[atom1->atomic_number];
            std::string element2 = atomicSymbols[atom2->atomic_number];

            if (is_bonded(*atom1,*atom2))
            {
                Bond tmp = {atom1, atom2, calc_distance(*atom1,*atom2), 0.0};
                bonds.push_back(tmp);
            }
        }
    }
    return bonds;
}

std::vector<Angle> GetAnglesFromBonds(std::vector<Atom> &atoms, std::vector<Bond> &bonds)
{
    std::vector<Angle> angles;
    for (unsigned int i = 0; i < bonds.size()-1; i++)
    {
        for (unsigned int j = i+1; j < bonds.size(); j++)
        {
            Atom* b1a1 = bonds[i].atom_1;
            Atom* b1a2 = bonds[i].atom_2;
            Atom* b2a1 = bonds[j].atom_1;
            Atom* b2a2 = bonds[j].atom_2;
            int iA = bonds[i].atom_1->index;
            int iB = bonds[i].atom_2->index;
            int iC = bonds[j].atom_1->index;
            int iD = bonds[j].atom_2->index;
            if (iA == iC)
            {   
                if (iB < iD)
                {
                    Angle tmp = {b1a2,b1a1,b2a2,calc_angle(*b1a2,*b1a1,*b2a2)};
                    angles.push_back(tmp);
                }
                else
                {
                    Angle tmp = {b2a2,b1a1,b1a2,calc_angle(*b1a2,*b1a1,*b2a2)};
                    angles.push_back(tmp);
                }
            }
            else if (iA == iD)
            {
                if (iB < iC)
                {
                    Angle tmp = {b1a2,b1a1,b2a1,calc_angle(*b1a2,*b1a1,*b2a1)};
                    angles.push_back(tmp);
                }
                else
                {
                    Angle tmp = {b2a1,b1a1,b1a2,calc_angle(*b1a2,*b1a1,*b2a1)};
                    angles.push_back(tmp);
                }                
            }
            else if (iB == iC)
            {
                if (iA < iD)
                {
                    Angle tmp = {b1a1,b1a2,b2a2,calc_angle(*b1a1,*b1a2,*b2a2)};
                    angles.push_back(tmp);
                }
                else
                {
                    Angle tmp = {b2a2,b1a2,b1a1,calc_angle(*b1a1,*b1a2,*b2a2)};
                    angles.push_back(tmp);
                }                
            }
            else if (iB == iD)
            {
                if (iA < iC)
                {
                    Angle tmp = {b1a1,b1a2,b2a1,calc_angle(*b1a1,*b1a2,*b2a1)};
                    angles.push_back(tmp);
                }
                else
                {
                    Angle tmp = {b2a1,b1a2,b1a1,calc_angle(*b1a1,*b1a2,*b2a1)};
                    angles.push_back(tmp);
                }
            }
            else
            {
                continue;
            }
        }
    }
    return angles;
}

std::vector<Dihedral> GetDihedralsFromAngles(std::vector<Atom> &atoms, std::vector<Bond> &bonds, std::vector<Angle> &angles)
{
    std::vector<Dihedral> dihedrals;
    for (Angle ang : angles)
    {
        Atom* a1 = ang.atom_1;
        Atom* a2 = ang.atom_2;
        Atom* a3 = ang.atom_3;
        for (Bond bond : bonds)
        {
            Atom* b1 = bond.atom_1;
            Atom* b2 = bond.atom_2;
            // angle center
            bool b1a2=(b1->index == a2->index);
            bool b2a2=(b2->index == a2->index);
            if ((!b1a2) && (!b2a2)) // If neither atom in the bond is the center of the angle, move on immediately.
            {
                continue;
            }

            // angle ends
            
            // if bond atom 1 is the same as angle center, ba2 must NOT be either a1 or a3
            if (b1a2)
            {
                bool b2a1=(b2->index == a1->index);
                bool b2a3=(b2->index == a3->index);
                if ((b2a1) || (b2a3)) // if either a1 or a3 is the same as b2, bail.
                {
                    continue;
                }
                Dihedral tmp = {a1,a2,a3,b2,calc_torsion(*a1,*a2,*a3,*b2)};
                dihedrals.push_back(tmp);
            }
            if (b2a2)
            {
                bool b1a1=(b1->index == a1->index);
                bool b1a3=(b1->index == a3->index);
                if ((b1a1) || (b1a3)) // if either a1 or a3 is the same as b2, bail.
                {
                    continue;
                }
                Dihedral tmp = {a1,a2,a3,b1,calc_torsion(*a1,*a2,*a3,*b1)};
                dihedrals.push_back(tmp);
            }
        }
    }
    return dihedrals;
}

std::vector<Torsion> GetTorsionsFromAngles(std::vector<Atom> &atoms, std::vector<Bond> &bonds, std::vector<Angle> &angles)
{
    std::vector<Torsion> torsions;
    for (Angle ang : angles)
    {
        Atom* a1 = ang.atom_1;
        Atom* a2 = ang.atom_2;
        Atom* a3 = ang.atom_3;
        for (Bond bond : bonds)
        {
            Atom* b1 = bond.atom_1;
            Atom* b2 = bond.atom_2;
            // angle center
            bool b1a2=(b1->index == a2->index);
            bool b2a2=(b2->index == a2->index);
            if ((b1a2) || (b2a2)) // If either atom in the bond is the center of the angle, move on immediately.
            {
                continue;
            }

            // angle ends
            bool b2a1=(b2->index == a1->index);
            bool b2a3=(b2->index == a3->index);
            bool b1a1=(b1->index == a1->index);
            bool b1a3=(b1->index == a3->index);

            // if bond atom 1 is the same as angle center, ba2 must NOT be either a1 or a3
            if (b2a1)
            {
                Torsion tmp = {b1,a1,a2,a3,calc_torsion(*b1,*a1,*a2,*a3)};
                torsions.push_back(tmp);
            }
            else if (b2a3)
            {
                Torsion tmp = {a1,a2,a3,b1,calc_torsion(*a1,*a2,*a3,*b1)};
                torsions.push_back(tmp);
            }
            else if (b1a1)
            {
                Torsion tmp = {b2,a1,a2,a3,calc_torsion(*b2,*a1,*a2,*a3)};
                torsions.push_back(tmp);
            }
            else if (b1a3)
            {
                Torsion tmp = {a1,a2,a3,b2,calc_torsion(*a1,*a2,*a3,*b2)};
                torsions.push_back(tmp);
            }
        }
    }
    return torsions;
}

std::map<int,int> GetNumBondsPerAtom(std::vector<Bond> bonds)
{
    std::map<int,int> n_bonds_per_atom={};
    for (Bond bond : bonds)
    {
        Atom* A = bond.atom_1;
        Atom* B = bond.atom_2;
        if (n_bonds_per_atom.count(A->index) == 0)
        {
            n_bonds_per_atom[A->index] = 0;
        }
        n_bonds_per_atom[A->index]++;
        if (n_bonds_per_atom.count(B->index) == 0)
        {
            n_bonds_per_atom[B->index] = 0;
        }
        n_bonds_per_atom[B->index]++;
    }
    return n_bonds_per_atom;
}

bool is_halogen(Atom A)
{
    std::vector<int> halogens = {9,17,35,53,85,117};
    if (std::find(halogens.begin(),halogens.end(),A.atomic_number) != halogens.end())
    {
        return true;
    }
    return false;
}

bool is_ring_atom(Atom A, std::vector<Ring> rings)
{
    int idx = A.index;
    for (Ring ring : rings)
    {
        for (Atom* at : ring.atoms)
        {
            if (at->index == idx)
            {
                return true;
            }
        }
    }
    return false;
}

std::vector<Ring> GetRingsFromAll(std::vector<Atom> &atoms, std::vector<Bond> &bonds, std::vector<Angle> &angles, std::vector<Torsion> &torsions)
{
    std::vector<Ring> rings = {};
    
    for (Torsion tors : torsions)
    {
        Atom* tA = tors.atom_1;
        Atom* tB = tors.atom_2;
        Atom* tC = tors.atom_3;
        Atom* tD = tors.atom_4;

        // 4 member rings:  torsions with matching endpoints as bond pairs.  Reminder that bond pairs are given in numerical order of atom index.
        for (Bond bond : bonds)
        {   
            Atom* bA = bond.atom_1;
            Atom* bB = bond.atom_2;
            if ((bA->index == tA->index) && (bB->index == tD->index))
            {
                Ring tmp={{tA,tB,tC,tD}};
                rings.push_back(tmp);
                continue;
            }
            if ((bB->index == tA->index) && (bA->index == tD->index))
            {
                Ring tmp={{tD,tC,tB,tA}};
                rings.push_back(tmp);
                continue;
            }
        }
        
        // 5 member rings:  torsions with matching endpoints in angles.
        for (Angle ang : angles)
        {
            Atom* aA = ang.atom_1;
            Atom* aB = ang.atom_2;
            Atom* aC = ang.atom_3;
            if ((aA->index == tA->index) && (aC->index == tD->index))
            {
                Ring tmp={{tA,tB,tC,tD,aB}};
                rings.push_back(tmp);
                continue;
            }
            if ((aC->index == tA->index) && (aA->index == tD->index))
            {
                Ring tmp={{tD,tC,tB,tA,aB}};
                rings.push_back(tmp);
                continue;
            }
        }
        // 6 member rings:  torsions with matching other torsions.
        for (Torsion ators : torsions)
        {
            Atom* atA = ators.atom_1;
            Atom* atB = ators.atom_2;
            Atom* atC = ators.atom_3;
            Atom* atD = ators.atom_4;
            if (((atA->index == tA->index) && (atD->index == tD->index)) || ((atA->index == tD->index) && (atD->index == tA->index)))
            { // if the endpoints match, we might have a ring.
                if (((atB->index == tB->index) && (atC->index == tC->index)) || ((atC->index == tB->index) && (atB->index == tC->index)))
                {   // if the MIDDLE points also match, we're just looking at the same torsion twice.
                    continue;
                }
                Ring tmp = {{tD,tC,tB,tA,atB,atC}};
                rings.push_back(tmp);
            }
        }
    }
    return rings;
}

void ObtainBondOrders(std::vector<Atom> &atoms, std::vector<Bond> &bonds, std::vector<Angle> angles)
{
    // Here's where we should be able to pass a collection of atoms, bonds, and angles, then determine the bond order for each.
    // The first stage should be running through all bonds with hydrogens - these will always have a bond order of 1 (single bonds).
    // Next is running through bonds containing carbons.  Any carbon with four total bonds will necessarily have all single bonds.
    // Since we are passing the bonds by reference, we can modify the bond orders directly (Bond struct has subunit 'double bond_order').
    
    // Get number of bonds for each atom.
    std::map<int,int> n_bonds_per_atom = GetNumBondsPerAtom(bonds);
    std::map<int,int> bond_orders_assigned = {};
    // Now check bond orders one-by-one.
    for (Bond &bond : bonds)
    {
        // if the bond order has already been assigned, move on.
        if (bond.bond_order > 0)
        {
            continue;
        }
        Atom* A = bond.atom_1;
        Atom* B = bond.atom_2;
        // Initialize bond_orders_assigned if not already present.
        if (bond_orders_assigned.count(A->index) == 0)
        {
            bond_orders_assigned[A->index] = 0;
        }
        if (bond_orders_assigned.count(B->index) == 0)
        {
            bond_orders_assigned[B->index] = 0;
        }

        // Check for hydrogens first.  This is generally a significant proportion of the total number of bonds, so doing this first will save runtime.
        if ((A->atomic_number == 1) || (B->atomic_number == 1))
        {
            bond.bond_order = 1.0;
            bond_orders_assigned[A->index]++;
            bond_orders_assigned[B->index]++;
            continue;
        }
        // Check for halogens.
        if ((is_halogen(*A)) || (is_halogen(*B)))
        {
            bond.bond_order = 1.0;
            bond_orders_assigned[A->index]++;
            bond_orders_assigned[B->index]++;
            continue;
        }
        // Check for sp3 carbons
        if (((A->atomic_number == 6) && (n_bonds_per_atom[A->index] == 4)) || ((B->atomic_number == 6) && (n_bonds_per_atom[B->index] == 4)))
        {
            bond.bond_order = 1.0;    
            bond_orders_assigned[A->index]++;
            bond_orders_assigned[B->index]++;
            continue;
        }
        // Check for sp3 nitrogens
        if (((A->atomic_number == 7) && (n_bonds_per_atom[A->index] == 4)) || ((B->atomic_number == 7) && (n_bonds_per_atom[B->index] == 4)))
        {
            bond.bond_order = 1.0;    
            bond_orders_assigned[A->index]++;
            bond_orders_assigned[B->index]++;
            continue;
        }

        // Check for oxygens next.
        if (A->atomic_number == 8)
        {
            // if 2 other atoms are bonded to this oxygen, then the bond must be single.
            if (n_bonds_per_atom[A->index] > 1)
            {
                bond.bond_order = 1.0;
                bond_orders_assigned[A->index]++;
                bond_orders_assigned[B->index]++;
                continue;
            }
            // If only one other atom is bonded, this can be a carbonyl (sp2 carbon) or a deprotonated hydroxyl group ()
            if ((B->atomic_number == 6) && (n_bonds_per_atom[B->index] == 4)) // bonded to sp3 carbon
            {
                bond.bond_order = 1.0;
                bond_orders_assigned[A->index]++;
                bond_orders_assigned[B->index]++;
                continue;
            }
        }
        else if (B->atomic_number == 8)
        {
            
        }
    }
}


// Molecule Class functions
//Overloaded Constructors
Molecule::Molecule()
{
    // empty Molecule, must add atoms manually from this point.
}

Molecule::Molecule(std::string filename)
{
    // identify filetype, then use internal parsers as appropriate.

}

Molecule::Molecule(std::string filename1,std::string filename2)
{
    // specifically useful if given two files (prmtop/rst7), identify which is which and then pass to specific parser.

}

Molecule::Molecule(std::vector<std::string> mol_chunk)
{
    // used when splitting up a multi-molecule PDB into separate pieces.

}

//Destructor
Molecule::~Molecule()
{

}

// Manually add atoms one-by-one or get called by constructor with arguments.
void Molecule::AddAtom(Atom& A)
{
    // check if Atom A already exists (check xyz values.  If there's a match, don't add atom and alert the user.)

    // Add atom to atoms vector.
    atoms.push_back(A);

    // Check for bonds by examining distances between all other atoms.
    for (Atom B : atoms)
    {
        if (is_bonded(A, B))
        {
            _add_bond(&A, &B);
        }
    }
}

void Molecule::_add_bond(Atom *A, Atom *B)
{
    int idx1 = A->index;
    int idx2 = B->index;
    for (Bond b : bonds) // make sure the bond doesn't already exist in this system.
    {
        int idx3 = b.atom_1->index;
        int idx4 = b.atom_2->index;
        if ((idx1 == idx3) && (idx2 == idx4))
        {
            return;
        }
    }

    if (idx1 < idx2) // sort from lowest index to highest, then add the bond.
    {
        double bondlength = calc_distance(*A,*B);
        Bond newbond = {A,B,bondlength,0.0};
        bonds.push_back(newbond);
    }
    else
    {
        double bondlength = calc_distance(*A,*B);
        Bond newbond = {A,B,bondlength,0.0};
        bonds.push_back(newbond);
    }

    // Check for any new angles formed by adding bond.
    for (Bond b : bonds)
    {
        int idx3 = b.atom_1->index;
        int idx4 = b.atom_2->index;
        if ((idx1 == idx3) || (idx1 == idx4)) // atom A is the vertex
        {
            if ((idx2 == idx4)||(idx2 == idx3))
            {
                continue; // Both indices match, so it's just the same bond again
            } 
            _add_angle(B,A,b.atom_2);
        }
        else if ((idx2 == idx3) || (idx2 == idx4)) // atom B is the vertex
        {
            if ((idx1 == idx4) || (idx1 == idx3))
            {
                continue; // Both indices match, so it's just the same bond again
            }
            _add_angle(A,B,b.atom_2);
        }
    }
}

void Molecule::_add_angle(Atom *A, Atom *B, Atom *C)
{
    int idx1 = A->index;
    int idx2 = B->index;
    int idx3 = C->index;

    for (Angle a : angles) // Check if angle is already in the molecule.
    {
        int idx_a = a.atom_1->index;
        int idx_b = a.atom_2->index;
        int idx_c = a.atom_3->index;
        if ((idx_a == idx1) && (idx_b == idx2) && (idx_c == idx3))
        {
            return;
        }
        if ((idx_c == idx1) && (idx_b == idx2) && (idx_a == idx3))
        {
            return;
        }
    }
    // get smaller endpoint index
    if (idx1 < idx3)
    {   
        double degrees = calc_angle(*A,*B,*C);
        Angle ang = {A,B,C,degrees};
        angles.push_back(ang);
    }
    else
    {   
        double degrees = calc_angle(*C,*B,*A);
        Angle ang = {C,B,A,degrees};
        angles.push_back(ang);
    }

    // check for formed dihedrals and torsions
    for (Bond b : bonds)
    {
        int idx4 = b.atom_1->index;
        int idx5 = b.atom_2->index;
        if ((angles.back().atom_1->index == idx4) && (angles.back().atom_2->index != idx5)) // torsion:  b2, a1, a2, a3
        {   
            _add_torsion(b.atom_2,angles.back().atom_1,angles.back().atom_2,angles.back().atom_3);
        }
        else if ((angles.back().atom_3->index == idx4) && (angles.back().atom_2->index != idx5))// torsion:  a1, a2, a3, b2
        {
            _add_torsion(angles.back().atom_1,angles.back().atom_2,angles.back().atom_3,b.atom_2);
        }
        else if ((angles.back().atom_1->index == idx5) && (angles.back().atom_2->index != idx4))// torsion:  b1, a1, a2, a3
        {
            _add_torsion(b.atom_1,angles.back().atom_1,angles.back().atom_2,angles.back().atom_3);
        }
        else if ((angles.back().atom_3->index == idx5) && (angles.back().atom_2->index != idx4))// torsion:  a1, a2, a3, b1
        {
            _add_torsion(angles.back().atom_1,angles.back().atom_2,angles.back().atom_3,b.atom_1);
        }
        else if ((angles.back().atom_2->index == idx4) && (angles.back().atom_1->index != idx5) && (angles.back().atom_3->index != idx5)) // dihedral: a1,a2,a3, b2
        {
            _add_dihedral(angles.back().atom_1,angles.back().atom_2,angles.back().atom_3,b.atom_2);
        }
        else if ((angles.back().atom_2->index == idx5) && (angles.back().atom_1->index != idx4) && (angles.back().atom_3->index != idx4)) // dihedral: a1,a2,a3, b1
        {
            _add_dihedral(angles.back().atom_1,angles.back().atom_2,angles.back().atom_3,b.atom_1);
        }
    }
}

void Molecule::_add_torsion(Atom* A, Atom* B, Atom* C, Atom* D)
{
    int idx_a = A->index;
    int idx_b = B->index;
    int idx_c = C->index;
    int idx_d = D->index;
    std::set <int> new_tors_indices = {idx_a,idx_b,idx_c,idx_d};
    for (Torsion tors : torsions)
    {
        int idx_1 = tors.atom_1->index;
        int idx_2 = tors.atom_2->index;
        int idx_3 = tors.atom_3->index;
        int idx_4 = tors.atom_4->index;
        std::set <int> tors_indices = {idx_1,idx_2,idx_3,idx_4};
        if (new_tors_indices == tors_indices)
        {
            return;
        }
    }

    if (idx_b < idx_c)
    {
        double degrees = calc_torsion(*A,*B,*C,*D);
        Torsion new_tors = {A,B,C,D,degrees};
        torsions.push_back(new_tors);
    }
    else
    {
        double degrees = calc_torsion(*D,*C,*B,*A);
        Torsion new_tors = {D,C,B,A,degrees};
        torsions.push_back(new_tors);
    }



}

void Molecule::_add_dihedral(Atom* A, Atom* B, Atom* C, Atom* D)
{
    // Assumes A-B-C is an angle and D is bonded to B as well.
    int idx_a = A->index;
    int idx_b = B->index;
    int idx_c = C->index;
    int idx_d = D->index;
    std::set <int> new_tors_indices = {idx_a,idx_b,idx_c,idx_d};
    for (Dihedral dihe : dihedrals)
    {
        int idx_1 = dihe.atom_1->index;
        int idx_2 = dihe.atom_2->index;
        int idx_3 = dihe.atom_3->index;
        int idx_4 = dihe.atom_4->index;
        std::set <int> tors_indices = {idx_1,idx_2,idx_3,idx_4};
        if (new_tors_indices == tors_indices)
        {
            return;
        }
    }
    double degrees = calc_dihedral(*A,*B,*C,*D);
    Dihedral dihe = {A,B,C,D,degrees};
    dihedrals.push_back(dihe);    
}


void Molecule::SavePDB(std::string filename)
{

}

void Molecule::SaveXYZ(std::string filename)
{
    std::stringstream buffer;
    buffer.str("");
    buffer << atoms.size() << std::endl; // XYZ format has the number of atoms first.
    buffer << " comment line " << std::endl; // XYZ format must have the comment line
    for (Atom atom : atoms)
    {
        // XYZ format has each atom on a single line, in the format
        // atomic_symbol   x_coordinate     y_coordinate    z_coordinate
        buffer << std::setw(4) << std::left << atomicSymbols[atom.atomic_number];
        buffer << std::setw(12) << std::right << std::setprecision(6) << atom.x;
        buffer << std::setw(12) << std::right << std::setprecision(6) << atom.y;
        buffer << std::setw(12) << std::right << std::setprecision(6) << atom.z;
        buffer << std::endl;
    }
    std::ofstream ofile(filename, std::ios::out);
    ofile << buffer.str();
    ofile.close();
}

void Molecule::SaveMOL2(std::string filename)
{

}
