#include "fileparsers.h"
#include "pdb.h"
#include "xyz.h"
#include "prmtop.h"
#include "mol2.h"

MolecularFile::MolecularFile(std::string filename)
{
    input_file = filename;
    filetype = fs::absolute(input_file).extension();
}

MolecularFile::~MolecularFile()
{
}

void MolecularFile::ParseFile()
{
    if (filetype == "pdb")
    {
        ParsePDB();
        return;
    }   
    if (filetype == "xyz")
    {
        ParseXYZ();
        return;
    }   
    if (filetype == "prmtop")
    {
        ParsePRMTOP();
        return;
    }   
    if (filetype == "mol2")
    {
        ParseMOL2();
        return;
    }
    error_log("Unable to parse given file.",1);
    return;
}
