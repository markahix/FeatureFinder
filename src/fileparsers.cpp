#include "fileparsers.h"
#include "pdb.h"

MolecularFile::MolecularFile(std::string filename)
{
    input_file = filename;
    filetype = utils::get_file_extension(input_file);
}

MolecularFile::~MolecularFile()
{
}

void MolecularFile::ParseFile()
{
    switch(filetype)
    {
        case 0:
            ParsePDB();
            break;
        case 1:
            ParseXYZ();
            break;
        case 2:
            ParsePRMTOP();
            break;
        case 3:
            ParseMOL2();
            break;
        default:
            std::cout << "Unable to parse given file." << std::endl;
            break;
    }
}

void MolecularFile::ParsePDB()
{
    PDB pdb(input_file);
    pdb.ReadPDB();
    
    // std::string line;
    // std::stringstream line_buf;
    // // std::ifstream pdb(input_file,std::ios::in);
    // std::vector<int> atoms_per_molecule = {};
    // std::vector<int> residues_per_molecule = {};
    // std::vector<std::string> resnames = {};
    // bool first_res_found = false;
    // int atom_count = 0;
    // int res_count = 0;
    // int old_res_num;
    // int new_res_num;
    // while (getline(pdb,line))
    // {
    //     if (line.find("END") != std::string::npos)
    //     {
    //         // Reached the end of the PDB's useful information, bail out.
    //         if (atom_count != 0)
    //         {
    //             atoms_per_molecule.push_back(atom_count);
    //             residues_per_molecule.push_back(res_count);
    //         }
    //         break;
    //     }
    //     if (line.find("TER") != std::string::npos)
    //     {
    //         atoms_per_molecule.push_back(atom_count);
    //         atom_count = 0;
    //         residues_per_molecule.push_back(res_count);
    //         res_count = 0;
    //         first_res_found = false;
    //     }
    //     if ((line.find("ATOM  ") != std::string::npos) || (line.find("HETATM") != std::string::npos))
    //     {
    //         atom_count++;
    //         if (!first_res_found)
    //         {
    //             old_res_num = stoi(line.substr(22,4));
    //             std::string resname = line.substr(17,3);
    //             resnames.push_back(resname);
    //             first_res_found = true;
    //             res_count++;
    //         }
    //         new_res_num = stoi(line.substr(22,4));
    //         if (new_res_num != old_res_num)
    //         {
    //             res_count++;
    //             std::string resname = line.substr(17,3);
    //             resnames.push_back(resname);
    //             old_res_num = new_res_num;
    //         }
    //     }
    // }

    // std::cout << "Atoms per molecule: ";
    // for (std::vector<int>::iterator it = atoms_per_molecule.begin(); it < atoms_per_molecule.end(); it++)
    // {
    //     std::cout << *it << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "Residues per molecule: ";
    // for (std::vector<int>::iterator it = residues_per_molecule.begin(); it < residues_per_molecule.end(); it++)
    // {
    //     std::cout << *it << " ";
    // }
    // std::cout << std::endl;

    // std::vector<std::string> molecule_types = {};
    // std::vector<std::string> molecule_masks = {};


}

void MolecularFile::ParseXYZ()
{

}

void MolecularFile::ParsePRMTOP()
{

}

void MolecularFile::ParseMOL2()
{

}
