#ifndef JOBSETTINGS_H
#define JOBSETTINGS_H

#include <string>
#include "utilities.h"

class JobSettings
{
    private:
        void ParseArgs(int argc, char **argv);
    public:
        // Class Variables
        fs::path input_file_name;
        bool debug_output;
        bool from_AGIMUS;
        fs::path AGIMUS_tracker_path;

        // Class Functions
        JobSettings(int argc, char **argv);
        ~JobSettings();
        
};

#endif