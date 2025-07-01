#ifndef FILEPARSERS_H
#define FILEPARSERS_H

#include "utilities.h"

class MolecularFile
{
    private:
        std::string input_file;
        int filetype;
        void ParsePDB();
        void ParseXYZ();
        void ParsePRMTOP();
        void ParseMOL2();

    public:
        MolecularFile(std::string filename);
        ~MolecularFile();
        void ParseFile();
};

#endif