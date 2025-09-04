#include "chemistry.h"

const std::vector<std::string> atomicSymbols = { "NULL",
    "H",  "HE", "LI", "BE", "B",  "C",  "N",  "O",  "F",  "NE",
    "NA", "MG", "AL", "SI", "P",  "S",  "CL", "AR", "K",  "CA",
    "SC", "TI", "V",  "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
    "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y",  "ZR",
    "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN",
    "SB", "TE", "I",  "XE", "CS", "BA", "LA", "CE", "PR", "ND",
    "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB",
    "LU", "HF", "TA", "W",  "RE", "OS", "IR", "PT", "AU", "HG",
    "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH",
    "PA", "U",  "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM",
    "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS",
    "RG", "CN", "NH", "FL", "MC", "LV", "TS", "OG"
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

int get_atomic_number(std::string symbol)
{
    for (int i = 0; i < atomicSymbols.size(); i++)
    {
        if (symbol == atomicSymbols[i])
        {
            return i;
        }
    }
}

double get_vdw_radius(int atomic_number)
{
    return vanDerWaalsRadii[atomic_number];
}

double get_vdw_radius(std::string symbol)
{
    return get_vdw_radius(get_atomic_number(symbol));
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
    double vdw1 = get_vdw_radius(atomicSymbols[atom1.atomic_number]);
    double vdw2 = get_vdw_radius(atomicSymbols[atom2.atomic_number]);
    double max_dist = (vdw1 + vdw2)/2;
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
    for (int i = 0; i < vec_ab.size(); i++)
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
    for (int i = 0; i < atoms.size()-1; i++ )
    {
        for (int j = i+1; j < atoms.size(); j++ )
        {
            Atom* atom1 = &atoms[i];
            Atom* atom2 = &atoms[j];
            std::string element1 = atomicSymbols[atom1->atomic_number];
            std::string element2 = atomicSymbols[atom2->atomic_number];

            if (is_bonded(*atom1,*atom2))
            {
                Bond tmp = {atom1, atom2, calc_distance(*atom1,*atom2), 1.0};
                bonds.push_back(tmp);
            }
        }
        return bonds;
    }
}

std::vector<Angle> GetAnglesFromBonds(std::vector<Atom> &atoms, std::vector<Bond> &bonds)
{
    std::vector<Angle> angles;
    for (int i = 0; i < bonds.size()-1; i++)
    {
        for (int j = i+1; j < bonds.size(); j++)
        {
            Atom* b1a1 = bonds[i].atom_1;
            Atom* b1a2 = bonds[i].atom_2;
            Atom* b2a1 = bonds[j].atom_1;
            Atom* b2a2 = bonds[j].atom_2;
            if (b1a1->index == b2a1->index)
            {   
                Angle tmp = {b1a2,b1a1,b2a2,calc_angle(*b1a2,*b1a1,*b2a2)};
                angles.push_back(tmp);
            }
            else if (b1a1->index == b2a2->index)
            {
                Angle tmp = {b1a2,b1a1,b2a1,calc_angle(*b1a2,*b1a1,*b2a1)};
                angles.push_back(tmp);
            }
            else if (b1a2->index == b2a1->index)
            {
                Angle tmp = {b1a1,b1a2,b2a2,calc_angle(*b1a1,*b1a2,*b2a2)};
                angles.push_back(tmp);
            }
            else if (b1a2->index == b2a2->index)
            {
                Angle tmp = {b1a1,b1a2,b2a1,calc_angle(*b1a1,*b1a2,*b2a1)};
                angles.push_back(tmp);
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