#ifndef FILEPARSERS_H
#define FILEPARSERS_H

#include "utilities.h"

class MolecularFile
{
    private:
        std::string input_file;
        std::string filetype;

    public:
        MolecularFile(std::string filename);
        ~MolecularFile();
        void ParseFile();
};

#endif