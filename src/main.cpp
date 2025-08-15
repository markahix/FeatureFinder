#include "utilities.h"
#include "jobsettings.h"
#include "fileparsers.h"

// First, identify type of input file
// - PDB:  identify proteins, nucleics, substrates, etc.
// - XYZ:  identify molecules based on connectivities, then if able, identify known substructures including amino acids, nucleic acids, etc.
// - MOL2: can't really think of anything needing to be done here if we've already got the MOL2 file...
// - prmtop:  identify proteins, nucleics, substrates, etc.

int main(int argc, char **argv)
{
    // Initialize variables
    JobSettings settings(argc,argv);

    // Obtain filename from command line
    std::cout << "Current input file: " << settings.input_file_name << std::endl;

    // Validate file is openable for reading and not empty.
    if (!validate_file(settings.input_file_name))
    {
        error_log("Input file errors were encountered.  Terminating FeatureFinder.",1);
    }
    
    // identify file type (prmtop, pdb, xyz, etc.)   
    MolecularFile mol_file(settings.input_file_name);
    mol_file.ParseFile();

    // Obtain feature types to identify (substructures of small molecules vs. proteins/nucleics/etc.)

    // identify features of molecules/systems accordingly.

    // output results to log in useful format for other modules.

    return 0;
}