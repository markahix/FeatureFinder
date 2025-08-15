#ifndef PDB_H
#define PDB_H

#include "utilities.h"

void ParsePDB();

/*
PDB can be divided up into several pieces:
Molecules -> Residues -> Atoms
Molecules are defined by connectivity and "TER" lines, Residues by number, name, and chain ID, and atoms by atom index.

Atoms just need their XYZ coordinates, elements, atom names (Amber types, etc.).  They can be stored under the residue to which they belong.

Residues just need their name, number, and chain ID, plus contain the collection of individual atoms within.  We also need to keep track of HEAD and TAIL atoms if that applies.

Molecules just need their number (internally assigned, not in the PDB itself) and the ORDERED collection of residues within.

Assumptions:
    Molecules will not have multiple chain IDs.  The presence of multiple chain IDs necessitates multiple molecules.
    Residues will have unique Name-Number-ChainID combinations.  
    Residues will have unique atom names internally.

Options for parsing a PDB:
    1.  Single-pass:  Read through the PDB one time, adding molecules, residues, and atoms line by line, but also needing to potentially work backwards to modify previous values in the event one line will change how we handle a previous line (such as head/tail atoms between residues).  Pros:  Probably more computationally efficient with file I/O.  Cons:  Might require more complicated work in handling the data we get.
    2.  Multi-pass:  Multiple reads through the PDB, allowing us to first assemble the data structure (pass for molecules, then pass for residues in molecules, then pass for atoms in residues), then fill in the details as we go back through.  Pros:  Might be easier to code.  Allows us to see if chain IDs are given inconsistently  Cons:  Might take more computational overhead.  Neutral: Might be negligible either way.

Let's look at the datastructures first.
Atoms:
    // Bookkeeping identifier
    int number;
    // location in space
    double x;
    double y;
    double z;
    // unique atomname, however if we move this to the residue level, this structure can be used for XYZs as well...  Or we just ... ignore it?
    std::string atomname;
    // element symbol (possibly do an enum with the periodic table?)
    std::string element;

Residues:
    
    // identifiers for residue.
    int number;
    std::string resname;
    std::string chainID; //If no chain ID given, assign "A"
    std::string head_atom;
    std::string tail_atom;
    std::vector<Atoms> atoms;

Molecules:
    int number;
    // If we use "chainID.resnum(leading zeros?).resname" as a unique identifier, then the Molecule level can have this as a map instead of a vector, and then we need only store the head, tail, and atoms vector...
    std::map<std::string, Residues> residues; //

PDBFile:
    std::string filename;
    std::vector<Molecules> molecules;

If these are all treated as structs, it might be easier to handle...  Then everything is just an instance of a given thing and we can modify it all internally...
This is strictly for parsing the PDB file and identifying connectivities.

Should we assume/assert that the user is providing a properly formatted PDB file?  Or should we attempt to handle error cases commonly seen?


Checking known structures:
Residues:  Known residue names (proteins, nucleics) can be checked against templates and alerts can be issued if they don't match up (in case of custom residues with names reserved for standards).  Additionally, residues that are known to be titratable can be identified, as can residues with specific protonation states like HIS/HIP/HID/HIE, CYS/CYM/CYX, ASP/ASH, GLU/GLH, etc.

    1.  First pass:  Does the residue have the same number of atoms as the template?
    2.  Second pass:  Does the residue have the same set of unique atom names as the template?
    3.  Third pass:  Are the atoms connected to each other in the same way (distance checks for all atom pairs)?

For non-matching structures, attempt to see if another match would be better and rename residue accordingly?
For unknown structures, extract and assign an AutoParams job?

*/

struct Atom
{
    int number;
    double x;
    double y;
    double z;
    std::string element;
    std::string atom_name;
};

struct Residue
{
    std::string head_atom;
    std::string tail_atom;
    std::vector<Atom> atoms;
};

struct Molecule
{
    std::map<std::string, Residue> residues;
};

struct PDBFile
{
    std::string filename;
    std::vector<Molecule> molecules;
};

void parse_pdb_line(std::string line, Molecule &molecule)
{
/*
ATOM      1  N   HIS A   1      49.668  24.248  10.436  1.00 25.00           N
......  ATOM/HETATM flag(1-6)
      ..... atom number (7-11)
           .
            .... atom name (13-16)
                . alternate location indicator
                 ... residue name (18-20)
                    . 
                     . chain identifier (22)
                      .... residue sequence number (23-26)
                          . code for insertion of residues
                           ...
                              ........ x coords (31-38)
                                      ........ y coords (39-46)
                                              ........ z coords (47-54)
                                                      ...... occupancy (55-60)
                                                            ...... temperature factor (61-66)
                                                                  .......
                                                                         ... segment identifier (73-76)
                                                                            .. elemental symbol (77-78)
                                                                              .. formal charge (79-80)
*/
    std::string atom_number  = line.substr(6 , 5);
    std::string atom_name    = line.substr(12, 4);
    std::string residue_name = line.substr(17, 3);
    std::string chain_id     = line.substr(21, 1);
    std::string residue_num  = line.substr(22, 4);
    std::string x_coords     = line.substr(30, 8);
    std::string y_coords     = line.substr(38, 8);
    std::string z_coords     = line.substr(46, 8);
    std::string occupancy    = line.substr(54, 6);
    std::string temperature  = line.substr(60, 6);
    std::string segment_id   = line.substr(72, 4);
    std::string element      = line.substr(76, 2);
    std::string charge       = line.substr(78, 2);
    // Assume "HETATM" or "ATOM  " lines only.
    // Parse line according to standard PDB formatting
    // If values don't align, throw an error and quit.
    // lines MUST have:  Residue name, residue number, atom name, x, y, z, element.
    // Useful:  chainID as well.
    bool BAD_DATA = false;
    if (is_empty(residue_name))
    {
        BAD_DATA = true;
    }
    if (is_empty(residue_num))
    {
        BAD_DATA = true;
    }
    if (is_empty(atom_name))
    {
        BAD_DATA = true;
    }
    if (is_empty(x_coords))
    {
        BAD_DATA = true;
    }
    if (is_empty(y_coords))
    {
        BAD_DATA = true;
    }
    if (is_empty(z_coords))
    {
        BAD_DATA = true;
    }
    if (is_empty(element))
    {
        BAD_DATA = true;
    }
    if (BAD_DATA)
    {
        error_log("Unable to parse line: " + line, 1);
    }

    if (is_empty(chain_id))
    {
        chain_id = "A";
    }
    // If we use "chainID.resnum(leading zeros?).resname" as a unique identifier, then the Molecule level can have this as a map instead of a vector, and then we need only store the head, tail, and atoms vector...
    // Generate Residue Key
    std::stringstream buffer;
    buffer.str("");
    buffer << chain_id << "." << std::setw(4) << std::setfill('0') << stoi(residue_num) << "." << residue_name;
    std::string reskey = buffer.str();
    if (molecule.residues.count(reskey) == 0)
    {
        molecule.residues[reskey] = {};
    }
    Atom new_atom {.number=stoi(atom_number), .x=stof(x_coords), .y=stof(y_coords), .z=stof(z_coords), .element = element, .atom_name = atom_name};
    molecule.residues[reskey].atoms.push_back(new_atom);
}


class PDBAtom
{
    private:
        int index;
        std::string atom_name;
        std::string residue_name;
        std::string chain_id;
        int residue_number;
        double x;
        double y;
        double z;
        std::string element;
        int charge;
    public:
        PDBAtom(std::string line);
        ~PDBAtom();
        void set_index(int idx);
        void set_atom_name(std::string aname);
        void set_residue_name(std::string rname);
        void set_chain_id(std::string chid);
        void set_residue_number(int rnum);
        void set_xyz(double newx, double newy, double newz);
        void set_element(std::string new_element);
        void set_charge(int new_charge);
        std::string get_atom_name();
};

class PDBResidue
{
    private:
        std::vector<PDBAtom> atoms;
        int resid;
        std::string resname;
        std::string chain_id;

    public:
        void AddAtom(std::string line);
        void CheckForErrors();
        void set_resid(int new_resid);
        void set_resname(std::string new_resname);
        void set_chain_id(std::string new_chid);
        PDBResidue(int id, std::string name, std::string chain);
        ~PDBResidue();
};

class PDBMolecule
{
    private:
        std::vector<PDBResidue> residues;
        int molecule_number;
    public:
        void AddResidue(int id, std::string name, std::string chain);
        PDBMolecule(int number);
        ~PDBMolecule();
};

class PDB
{
    private:
        std::string pdb_file;
        std::vector <std::string> molecule_ranges;
        std::vector <PDBMolecule> molecules;
        
    public:
        void ReadPDB();

        PDB(std::string filename);
        ~PDB();
};

#endif