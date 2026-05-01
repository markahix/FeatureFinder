#include "jobsettings.h"
JobSettings::JobSettings(int argc, char** argv)
{
    input_file_name = "file.pdb";
    debug_output = false;
    from_AGIMUS = false;
    AGIMUS_tracker_path = "./.AGIMUStracker";
    ParseArgs(argc, argv);
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

        if (curr_flag == "--fromAGIMUS")
        {
            from_AGIMUS = true;
            AGIMUS_tracker_path = fs::absolute(argv[i+1]);
            i+=2;
            continue;
        }
        
        else
        {
            std::cout << "Unrecognized keyword argument: " << argv[i] << std::endl;
            i++;
        }
    }
}