#include "jobsettings.h"
JobSettings::JobSettings()
{
    input_file_name = "file.pdb";
    debug_output = false;
}

JobSettings::~JobSettings()
{

}

void JobSettings::ParseArgs(int argc, char** argv)
{
    int i = 1;
    while (i < argc)
    {
        std::string curr_flag = argv[i];
        if (curr_flag == "--inputfile" || curr_flag == "-i")
        {
            input_file_name = fs::absolute(argv[i+1]);
            i+=2;
            continue;
        }

        if (curr_flag == "--debug")
        {
            debug_output = true;
            i++;
            continue;
        }

        else
        {
            std::cout << "Unrecognized keyword argument: " << argv[i] << std::endl;
            i++;
        }
    }
}