#include "utilities.h"
#include "jobsettings.h"
#include "fileparsers.h"

namespace fs = std::experimental::filesystem;

int main(int argc, char **argv)
{
    // Initialize variables
    JobSettings settings;

    // Obtain filename from command line
    settings.ParseArgs(argc,argv);
    std::cout << "Current input file: " << settings.input_file_name << std::endl;

    // Validate file is openable for reading and not empty.
    utils::ExitOnFailure(utils::file_exists_and_can_be_opened_and_is_not_empty(settings.input_file_name), "Input file errors were encountered.  Terminating FeatureFinder.");

    // identify file type (prmtop, pdb, xyz, etc.)   
    MolecularFile mol_file(settings.input_file_name);
    mol_file.ParseFile();

    // Obtain feature types to identify (substructures of small molecules vs. proteins/nucleics/etc.)

    // identify features of molecules/systems accordingly.

    // output results to log in useful format for other modules.

    return 0;
}