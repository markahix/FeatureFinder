#ifndef FILEPARSERS_H
#define FILEPARSERS_H

#include "utilities.h"
#include "jobsettings.h"

class Molecule
{
    public:
        Molecule(std::vector<std::string> mol_atom_names, std::vector<std::pair<int,int>> mol_edges);
        ~Molecule();
    private:
        std::vector<std::string> atom_names;
        std::vector<std::pair<int,int>> edges;
};

class MolecularFile
{
    private:
        std::string input_file;
        std::string filetype;
        JobSettings* settings_ptr;

    public:
        MolecularFile(std::string filename, JobSettings &settings);
        ~MolecularFile();
        void ParseFile();
};

/*
PDB Parser --> Molecule Format
MOL2 Parser --> Molecule Format
PRMTOP Parser --> Molecule Format
XYZ Parser --> Molecule Format
*/
// Molecule PDBtoMolecule(std::string filename, std::vector<std::string> &mol_atom_names, std::vector<std::pair<int,int>> &mol_edges);
// Molecule MOL2toMolecule(std::string filename, std::vector<std::string> &mol_atom_names, std::vector<std::pair<int,int>> &mol_edges);
// Molecule PRMTOPtoMolecule(std::string filename, std::vector<std::string> &mol_atom_names, std::vector<std::pair<int,int>> &mol_edges);
// Molecule XYZtoMolecule(std::string filename, std::vector<std::string> &mol_atom_names, std::vector<std::pair<int,int>> &mol_edges);
#endif