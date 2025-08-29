#include "pdb.h"


// PDBAtom Functions
PDBAtom::PDBAtom(std::string line)
{
    index = atoi(line.substr(6,5).c_str());
    atom_name = line.substr(12,4);
    residue_name = line.substr(17,3);
    chain_id = line.substr(21,1);
    if (chain_id == " ")
    {
        chain_id = "A";
    }
    residue_number = atoi(line.substr(22,4).c_str());
    x = atof(line.substr(30,8).c_str());
    y = atof(line.substr(38,8).c_str());
    z = atof(line.substr(46,8).c_str());
    element = line.substr(76,2);
    charge = 0;
    if (line.substr(78,2).find_first_not_of(' ') != std::string::npos)
    {
        charge = atoi(line.substr(78,2).c_str());
    }
}
PDBAtom::~PDBAtom()
{

}
void PDBAtom::set_index(int idx)
{
    index = idx;
}
void PDBAtom::set_atom_name(std::string aname)
{
    atom_name = aname;
}
void PDBAtom::set_residue_name(std::string rname)
{
    residue_name = rname;
}
void PDBAtom::set_chain_id(std::string chid)
{
    chain_id = chid;
}
void PDBAtom::set_residue_number(int rnum)
{
    residue_number = rnum;
}
void PDBAtom::set_xyz(double newx, double newy, double newz)
{
    x = newx;
    y = newy;
    z = newz;
}
void PDBAtom::set_element(std::string new_element)
{
    element = new_element;
}
void PDBAtom::set_charge(int new_charge)
{
    charge = new_charge;
}
std::string PDBAtom::get_atom_name()
{
    return atom_name;
}
//PDBResidue Functions

void PDBResidue::AddAtom(std::string line)
{
    atoms.push_back(PDBAtom(line));
}
void PDBResidue::CheckForErrors()
{   
    // Check for unique atom names.
    std::vector<std::string> tmp_atom_names = {};
    for (unsigned int i = 0;i < atoms.size(); i++)
    {
        tmp_atom_names.push_back(atoms[i].get_atom_name());
    }
    std::set<std::string> unique_names(tmp_atom_names.begin(), tmp_atom_names.end());
    if (unique_names.size() != tmp_atom_names.size())
    {
        std::cout << "ERROR: Residue " << resid << " may have non-unique atom names." << std::endl;
    }    
}

void PDBResidue::set_resid(int new_resid)
{
    resid = new_resid;
    for (unsigned int i = 0;i < atoms.size(); i++)
    {
        atoms[i].set_residue_number(resid);
    }
}
void PDBResidue::set_resname(std::string new_resname)
{
    resname = new_resname;
    for (unsigned int i = 0;i < atoms.size(); i++)
    {
        atoms[i].set_residue_name(resname);
    }
}
void PDBResidue::set_chain_id(std::string new_chid)
{
    chain_id = new_chid;
    for (unsigned int i = 0;i < atoms.size(); i++)
    {
        atoms[i].set_chain_id(chain_id);
    }
}

PDBResidue::PDBResidue(int id, std::string name, std::string chain)
{
    atoms = {};
    resid = id;
    resname = name;
    chain_id = chain;
}
PDBResidue::~PDBResidue()
{

}

// CONSTRUCTOR
PDB::PDB(std::string filename)
{
    pdb_file = filename;
    molecule_ranges = {};
}

// DESTRUCTOR
PDB::~PDB()
{
}

void PDB::ReadPDB()
{
    std::string line;
    std::ifstream ifile(pdb_file,std::ios::in);
    if (!ifile.is_open())
    {
        std::cout << "Unable to open " << pdb_file << ". Aborting. " << std::endl;
        return;
    }

    std::stringstream current_molecule_range;
    int current_residue = 0;
    int start_molecule = 0;
    bool first_atom_found = false;

    while (getline(ifile,line))
    {   
        if (line.find("END") != std::string::npos)
        {   // Reached the end of the PDB's useful information, bail out.
            break;
        }
        if (line.find("TER") != std::string::npos)
        {
            // End of molecule, close the current molecule and prepare another molecule.
            current_molecule_range.str("");
            current_molecule_range << ":" << start_molecule;
            if (current_residue != start_molecule)
            {
                current_molecule_range << "-" << current_residue;
            }
            molecule_ranges.push_back(current_molecule_range.str());
            first_atom_found = false;
            
            continue;
        }
        if ((line.find("ATOM") == std::string::npos) && (line.find("HETATM") == std::string::npos) )
        {   // If line is not an atom-containing line, continue.
            continue;
        }
        current_residue = atoi(line.substr(22,4).c_str());
        if (!first_atom_found)
        {
            start_molecule = atoi(line.substr(22,4).c_str());
            first_atom_found = true;
        }
    }
}
PDBMolecule::PDBMolecule(int number)
{

}

PDBMolecule::~PDBMolecule()
{

}

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
//OneLetterAminoAcids
std::map<std::string,std::string> OLAA =   {{"ALA","A"},{"CYS","C"},{"CYM","C"},{"CYX","C"},{"ASP","D"},{"ASH","D"},{"GLU","E"},
                                                {"GLH","E"},{"PHE","F"},{"GLY","G"},{"HIS","H"},{"HID","H"},{"HIE","H"},{"HIP","H"},
                                                {"ILE","I"},{"LYS","K"},{"LEU","L"},{"MET","M"},{"ASN","N"},{"PRO","P"},{"GLN","Q"},
                                                {"ARG","R"},{"SER","S"},{"THR","T"},{"VAL","V"},{"TRP","W"},{"TYR","Y"}};

std::map<std::string,int> AAcharges =   {{"ALA",0},{"CYS",0},{"CYM",-1},{"CYX",0},{"ASP",-1},{"ASH",0},{"GLU",-1},
                                                {"GLH",0},{"PHE",0},{"GLY",0},{"HIS",0},{"HID",0},{"HIE",0},{"HIP",1},
                                                {"ILE",0},{"LYS",1},{"LEU",0},{"MET",0},{"ASN",0},{"PRO",0},{"GLN",0},
                                                {"ARG",1},{"SER",0},{"THR",0},{"VAL",0},{"TRP",0},{"TYR",0}};


                                                //OneLetterNucleicAcids
std::map<std::string,std::string> OLNA =   {{"DA3","A"},{"DA ","A"},{"DA5","A"},{"A3 ","A"},{"A  ","A"},{"A5 ","A"},
                                            {"DC3","C"},{"DC ","C"},{"DC5","C"},{"C3 ","C"},{"C  ","C"},{"C5 ","C"},
                                            {"DG3","G"},{"DG ","G"},{"DG5","G"},{"G3 ","G"},{"G  ","G"},{"G5 ","G"},
                                            {"DT3","T"},{"DT ","T"},{"DT5","T"},{"U3 ","U"},{"U  ","U"},{"U5 ","U"}};

// Metal Ions List
std::vector<std::string> KnownMetalIons = {"ZN "};

// Counterions
std::vector<std::string> KnownCounterions = {"ZN "};

// Known Solvents
std::vector<std::string> KnownSolvents = {"WAT"};

bool is_solvents(std::vector<std::string> residue_keys)
{
    int n_residues = residue_keys.size();
    int prot_residues = 0;
    double percentage;
    for (std::string reskey : residue_keys)
    {
        std::string resname = reskey.substr(reskey.find_last_of('.')+1,reskey.size() - reskey.find_last_of('.') - 1);
        if (std::find(KnownSolvents.begin(),KnownSolvents.end(),resname) == KnownSolvents.end())
        {
            continue;
        }
        prot_residues++;
    }
    percentage = (double)prot_residues / (double)n_residues;
    if (percentage > .5)
    {
        return true;
    }
    return false;
}

bool is_counterions(std::vector<std::string> residue_keys)
{
    int n_residues = residue_keys.size();
    int prot_residues = 0;
    double percentage;
    for (std::string reskey : residue_keys)
    {
        std::string resname = reskey.substr(reskey.find_last_of('.')+1,reskey.size() - reskey.find_last_of('.') - 1);
        if (std::find(KnownCounterions.begin(),KnownCounterions.end(),resname) == KnownCounterions.end())
        {
            continue;
        }
        prot_residues++;
    }
    percentage = (double)prot_residues / (double)n_residues;
    if (percentage > .5)
    {
        return true;
    }
    return false;
}

bool is_metal_ions(std::vector<std::string> residue_keys)
{
    int n_residues = residue_keys.size();
    int prot_residues = 0;
    double percentage;
    for (std::string reskey : residue_keys)
    {
        std::string resname = reskey.substr(reskey.find_last_of('.')+1,reskey.size() - reskey.find_last_of('.') - 1);
        if (std::find(KnownMetalIons.begin(),KnownMetalIons.end(),resname) == KnownMetalIons.end())
        {
            continue;
        }
        prot_residues++;
    }
    percentage = (double)prot_residues / (double)n_residues;
    if (percentage > .5)
    {
        return true;
    }
    return false;

}

bool is_mostly_nucleic(std::vector<std::string> residue_keys)
{
    int n_residues = residue_keys.size();
    int prot_residues = 0;
    double percentage;
    for (std::string reskey : residue_keys)
    {
        std::string resname = reskey.substr(reskey.find_last_of('.')+1,reskey.size() - reskey.find_last_of('.') - 1);
        if (OLNA.count(resname) == 0)
        {
            continue;
        }
        prot_residues++;
    }
    percentage = (double)prot_residues / (double)n_residues;
    if (percentage > .5)
    {
        return true;
    }
    return false;
}

bool is_mostly_protein(std::vector<std::string> residue_keys)
{
    int n_residues = residue_keys.size();
    int prot_residues = 0;
    double percentage;
    for (std::string reskey : residue_keys)
    {
        if (OLAA.count(reskey.substr(reskey.find_last_of('.')+1,3)) == 0)
        {
            continue;
        }
        prot_residues++;
    }
    percentage = (double)prot_residues / (double)n_residues;
    if (percentage > .5)
    {
        return true;
    }
    return false;
}

std::string RawProteinSequence(std::vector<std::string> residue_keys)
{
    std::string sequence = "";
    for (std::string reskey : residue_keys)
    {
        std::string resname = reskey.substr(reskey.find_last_of('.')+1,reskey.size() - reskey.find_last_of('.') - 1);
        if (OLAA.count(resname) == 0)
        {
            OLAA[resname]="_";
        }
        sequence+= OLAA[resname];
    }
    return sequence;
}

int GetProteinTotalCharge(std::vector<std::string> residue_keys)
{
    int total_charge=0;
    for (std::string reskey : residue_keys)
    {
        std::string resname = reskey.substr(reskey.find_last_of('.')+1,reskey.size() - reskey.find_last_of('.') - 1);
        if (AAcharges.count(resname) == 0)
        {
            AAcharges[resname]=0;
            std::cout << "Unknown charge for " << resname << "." <<std::endl;
        }
        total_charge += AAcharges[resname];
    }
    return total_charge;
}

void Protein_Sequence_to_FASTA(std::vector<std::string> residue_keys, int molecule_number)
{
    int n_res = 1;
    int n_mol = 0;
    std::string prev_chain = residue_keys[0].substr(0,residue_keys[0].find_first_of('.'));
    std::stringstream fasta;
    fasta.str("");
    fasta << "Protein sequence detected." << std::endl;
    fasta << "Chain ID: " << prev_chain << std::endl;
    fasta << std::setw(4) << std::right << n_res << "  ";
    for (std::string reskey : residue_keys)
    {
        std::string curr_chain = reskey.substr(0,reskey.find_first_of('.'));
        if (curr_chain != prev_chain)
        {
            fasta << "\n\n";
            n_res = 1;
            n_mol++;
            fasta << "Chain ID: " << curr_chain << std::endl;
            fasta << std::setw(4) << std::right << n_res << "  ";
            prev_chain = curr_chain;
        }
        std::string resname = reskey.substr(reskey.find_last_of('.')+1,3);
        if (OLAA.count(resname) == 0)
        {
            OLAA[resname]="_";
        }
        fasta << OLAA[resname];
        if (n_res%10 == 0)
        {
            fasta << " ";
        }
        if (n_res %50 == 0)
        {
            fasta << " " << std::setw(4) << std::left << n_res << "\n";
            fasta << std::setw(4) << std::right << n_res+1 << "  ";
        }
        n_res++;
    }
    std::stringstream buffer;
    buffer.str("");
    buffer << "Protein_"<< std::setw(4) << std::right << std::setfill('0') << molecule_number << ".fasta";
    std::ofstream ofile;
    ofile.open(buffer.str(),std::ios::out);
    ofile << fasta.str() << std::endl;
    ofile.close();
}

void Nucleic_Sequence_to_FASTA(std::vector<std::string> residue_keys, int molecule_number)
{
    int n_res = 1;
    int n_mol = 0;
    std::string prev_chain = residue_keys[0].substr(0,residue_keys[0].find_first_of('.'));
    std::stringstream fasta;
    fasta.str("");
    fasta << "Nucleic acid sequence detected." << std::endl;
    fasta << "Chain ID: " << prev_chain << std::endl;
    fasta << std::setw(4) << std::right << n_res << "  ";
    for (std::string reskey : residue_keys)
    {
        std::string curr_chain = reskey.substr(0,reskey.find_first_of('.'));
        if (curr_chain != prev_chain)
        {
            fasta << "\n\n";
            n_res = 1;
            n_mol++;
            fasta << "Chain ID: " << curr_chain << std::endl;
            fasta << std::setw(4) << std::right << n_res << "  ";
            prev_chain = curr_chain;
        }
        std::string resname = reskey.substr(reskey.find_last_of('.')+1,reskey.size() - reskey.find_last_of('.') - 1);
        if (OLNA.count(resname) == 0)
        {
            OLNA[resname]="_";
        }
        fasta << OLNA[resname];
        if (n_res%10 == 0)
        {
            fasta << " ";
        }
        if (n_res %50 == 0)
        {
            fasta << " " << std::setw(4) << std::left << n_res << "\n";
            fasta << std::setw(4) << std::right << n_res+1 << "  ";
        }
        n_res++;
    }
    std::stringstream buffer;
    buffer.str("");
    buffer << "Nucleic_"<< std::setw(4) << std::right << std::setfill('0') << molecule_number << ".fasta";
    std::ofstream ofile;
    ofile.open(buffer.str(),std::ios::out);
    ofile << fasta.str() << std::endl;
    ofile.close();
}

std::string reskey_from_pdbline(std::string pdbline, std::string default_chain_id)
{
    std::string residue_name = pdbline.substr(17, 3);
    std::string chain_id     = pdbline.substr(21, 1);
    std::string residue_num  = pdbline.substr(22, 4);
    if (is_empty(chain_id))
    {
        chain_id = default_chain_id;
    }
    std::stringstream buffer;
    buffer.str("");
    buffer << chain_id << "." << std::setw(4) << std::setfill('0') << stoi(residue_num) << "." << residue_name;
    return buffer.str();
}

std::vector<std::vector<std::string>> PDBFile_to_Molecule_Sections(std::string inputfile)
{
    std::ifstream ifile(inputfile,std::ios::in);
    if (! ifile.is_open())
    {
        error_log("Unable to open input file: " + inputfile, 1);
    }

    std::string line;
    std::vector<std::string> residue_keys = {};

    std::string def_chain_id = "A";
    std::string prev_chain_id = "";

    std::vector<std::vector<std::string>> molecules = {{}}; // MAIN DATA TO BE RETURNED

    while(getline(ifile,line))
    {
        bool end_of_molecule = false;
        bool end_of_file     = false;
        bool skip_parse      = false;
        // If line has "END", it's the end of the molecule and end of the file.
        if (line.find("END") != std::string::npos)
        {
            end_of_molecule = true;
            end_of_file     = true;
            skip_parse      = true;
        }

        else if (line.find("TER") != std::string::npos)
        {
            end_of_molecule = true;
            skip_parse      = true;
        }

        else if ((line.find("HETATM") == std::string::npos) && (line.find("ATOM  ") == std::string::npos))
        {   // If the line doesn't have HETATM or ATOM in it, assume it's not information we need right now and move on.
            skip_parse      = true;
        }

        if (!skip_parse)
        {   
            std::string curr_chain_id = line.substr(21, 1);
            if (prev_chain_id == "")
            {
                prev_chain_id = curr_chain_id;
            }
            if (curr_chain_id != prev_chain_id)
            {
                end_of_molecule = true;
                prev_chain_id = curr_chain_id;
            }
        }
        
        if (end_of_molecule)
        {
            if (residue_keys.size() == 0)
            {
                continue;
            }
            residue_keys = {};
            molecules.push_back({});
        }

        if (!skip_parse)
        {
            molecules.back().push_back(line);
            std::string reskey = reskey_from_pdbline(line, def_chain_id);
            if (std::find(residue_keys.begin(), residue_keys.end(), reskey) == residue_keys.end())
            {
                residue_keys.push_back(reskey);
            }
        }

        if (end_of_file)
        {
            break;
        }
    }
    ifile.close();
    return molecules;
}

std::vector<std::string> GetReskeysFromMolecule(std::vector<std::string> mol_chunk, std::string def_chain_id)
{
    std::vector<std::string> residue_keys = {};
    for (std::string line : mol_chunk)
    {
        std::string reskey = reskey_from_pdbline(line, def_chain_id);
        if (std::find(residue_keys.begin(), residue_keys.end(), reskey) == residue_keys.end())
        {
            residue_keys.push_back(reskey);
        }
    }
    return residue_keys;
}

void SequenceMolecule(std::vector<std::string> residue_keys, int molecule_number)
{
        if (residue_keys.size() == 0)
        {
            return;
        }
        if (is_mostly_protein(residue_keys))
        {
            Protein_Sequence_to_FASTA(residue_keys, molecule_number);
        }
        else if (is_mostly_nucleic(residue_keys))
        {
            Nucleic_Sequence_to_FASTA(residue_keys,  molecule_number);
        }
        else if (is_metal_ions(residue_keys))
        {
            std::cout << "Metal Ions detected in molecule " <<  molecule_number << std::endl;
        }
}

void ParsePDB(std::string input_file)
{
    std::string def_chain_id = "A";
    std::vector<std::vector<std::string>> molecule_chunks = PDBFile_to_Molecule_Sections(input_file);
    // Handling for Multiple Molecules
    if (molecule_chunks.size() > 1)
    {
        for (unsigned int i = 0; i < molecule_chunks.size(); i ++)
        {
            std::vector<std::string> residue_keys = GetReskeysFromMolecule(molecule_chunks[i], def_chain_id);
            SequenceMolecule(residue_keys, i+1);
        }
    }
    // Handling for single-molecule
    else
    {
        // Get list of residues from the single molecule.
        // If the list of residues is just one residue, treat this as a small-molecule feature examination.
        // otherwise, sequence as normal.

    }
}

