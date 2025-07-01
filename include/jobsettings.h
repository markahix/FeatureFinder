#ifndef JOBSETTINGS_H
#define JOBSETTINGS_H

#include <string>
#include "utilities.h"

class JobSettings
{
    public:
        // Class Variables
        fs::path input_file_name;
        bool debug_output;

        // Class Functions
        JobSettings();
        ~JobSettings();
        void ParseArgs(int argc, char **argv);
};

#endif