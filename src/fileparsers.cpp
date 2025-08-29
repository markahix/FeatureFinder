#include "fileparsers.h"
#include "pdb.h"
#include "xyz.h"
#include "prmtop.h"
#include "mol2.h"

MolecularFile::MolecularFile(std::string filename)
{
    input_file = filename;
    filetype = fs::absolute(input_file).extension();
    std::cout << filetype << std::endl;
}

MolecularFile::~MolecularFile()
{
}

void MolecularFile::ParseFile()
{
    if (filetype == ".pdb")
    {
        normal_log("Given input file is a PDB.");
        ParsePDB(input_file);
        return;
    }   
    // if (filetype == "xyz")
    // {
    //     normal_log("Given input file is an XYZ.");
    //     ParseXYZ(input_file);
    //     return;
    // }   
    // if (filetype == "prmtop")
    // {
    //     ParsePRMTOP(input_file);
    //     normal_log("Given input file is a PRMTOP.");
    //     return;
    // }   
    // if (filetype == "mol2")
    // {
    //     normal_log("Given input file is an MOL2.");
    //     ParseMOL2(input_file);
    //     return;
    // }
    error_log("Unable to parse given file.",1);
    return;
}
