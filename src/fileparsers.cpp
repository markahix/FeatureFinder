#include "fileparsers.h"
#include "pdb.h"
#include "xyz.h"
#include "prmtop.h"
#include "mol2.h"

/*
TO DO:
Build centralized graph structure to hold atoms, bonds, and molecules (groups) from any/all accepted filetypes.

Atoms will need to have x,y,z, and element information at minimum.
// std::vector<Atom*> atoms;

Bonds will simply contain two Atom* instances (and a distance? or bond order?)
// std::vector<Bond*> bonds;

Molecules will have lists of atom and bond pointers.  Do I need an identifier?
// std::vector<Molecule*> molecules;

For multi-molecule runs, FeatureFinder will simply be identifying the individual molecules as we currently have them.

For single-molecule runs, FeatureFinder will identify specific features like functional groups or known scaffolds

I need a way to compare graphs of elements to see if a user-provided graph matches any known subgraphs

If I can get a graph comparison protocol established, then I can just check an incoming molecule against the list of known (and coded) subgraphs.
1. take subgraph as set of nodes.
2. see if all nodes in subgraph have matches in big graph - compare node weights (element) and number of edges with that node (# bonds)
    a. if not, big graph doesn't contain small graph, move on.
3. see if all edges in subgraph have matches in big graph
    a. if not, big graph doesn't contain small graph, move on.
4. keep comparing connections.  for each node in the subgraph, does the corresponding node in the big graph remain connected to the sample?
*/


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
    // if (filetype == ".xyz")
    // {
    //     normal_log("Given input file is an XYZ.");
    //     ParseXYZ(input_file);
    //     return;
    // }   
    // if (filetype == ".prmtop")
    // {
    //     ParsePRMTOP(input_file);
    //     normal_log("Given input file is a PRMTOP.");
    //     return;
    // }   
    // if (filetype == ".mol2")
    // {
    //     normal_log("Given input file is a MOL2.");
    //     ParseMOL2(input_file);
    //     return;
    // }
    error_log("Unable to parse given file.",1);
    return;
}
