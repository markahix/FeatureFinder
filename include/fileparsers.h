#ifndef FILEPARSERS_H
#define FILEPARSERS_H

#include "utilities.h"
#include "jobsettings.h"

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

#endif