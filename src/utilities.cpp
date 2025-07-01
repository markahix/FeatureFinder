
#include "utilities.h"
// #include "classes.h"
// #include "slurm.h"

namespace utils
{
    void splash_screen()
    {
        std::cout << "#############################" << std::endl;
        std::cout << "# Behold the Mighty AGIMUS! #" << std::endl;
        std::cout << "#  Deploying FeatureFinder  #" << std::endl;
        std::cout << "#############################" << std::endl;
    }

    void silent_shell(const char* cmd)
    {
        std::array<char, 128> buffer;
        std::string result;
        std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
        if (!pipe) throw std::runtime_error("popen() failed!");
        while (!feof(pipe.get())) {
            if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
                result += buffer.data();
        }
    }

    std::string GetSysResponse(const char* cmd)
    {
        std::array<char, 128> buffer;
        std::string result;
        std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
        if (!pipe) throw std::runtime_error("popen() failed!");
        while (!feof(pipe.get())) {
            if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
                result += buffer.data();
        }
        return result;
    }

    bool CheckProgAvailable(const char* program)
    {
            std::string result;
            std::string cmd;
            cmd = "which ";
            cmd += program;
            result=GetSysResponse(cmd.c_str());
            if (result.empty())
            {
                std::cout << "Missing program: " << program << std::endl;
                return false;
            }
            return true;
    }

    void write_to_file(std::string inputfilename, std::string buffer)
    {
        std::ofstream outFile;
        outFile.open(inputfilename,std::ios::out);
        outFile << buffer;
        outFile.close();
    }

    void append_to_file(std::string inputfilename, std::string buffer)
    {
        if (!CheckFileExists(inputfilename))
        {
            write_to_file(inputfilename,"");
        }
        std::ofstream outFile;
        outFile.open(inputfilename,std::ios::app);
        outFile << buffer;
        outFile.close();
    }


    bool CheckFileExists(std::string filename)
    {
        if ( fs::exists(filename) )
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
 
    std::string GetTimeAndDate()
    {
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::stringstream buffer;
        buffer.str("");
        buffer << std::put_time(&tm, "%Y.%m.%d %H:%M:%S");
        return buffer.str();
    }
 
    int is_empty(const char *s) 
    {
        while (*s != '\0') 
        {
            if (!isspace((unsigned char)*s))
            return 0;
            s++;
        }
        return 1;
    }
 
    std::string LastLineOfFile(std::string filename)
    {
        std::ifstream fin;
        fin.open(filename);
        std::string line;
        std::vector<std::string> lines_in_order;
        while (std::getline(fin, line))
        {   
            if ( !is_empty(line.c_str()) )
            {
                lines_in_order.push_back(line);
            }
        }
        fin.close();
        return lines_in_order[lines_in_order.size()-1];
    }

    int count_lines_in_file(std::string filename)
    {
        std::string dummy;
        int count = 0;
        std::ifstream file(filename);
        while (getline(file,dummy))
        {
            count ++;
        }
        return count;
    }

    std::string string_between(std::string incoming, std::string first_delim, std::string second_delim)
    {
        unsigned first = incoming.find(first_delim);
        if (first == std::string::npos)
        {
            return incoming;
        }
        unsigned last = incoming.find(second_delim);
        if (last == std::string::npos)
        {
            return incoming.substr(first + 1, incoming.size());
        }
        return incoming.substr(first + 1, last - first - 1);
    }
    
    void compress_and_delete(std::string directory)
    {
        std::stringstream buffer;
        buffer.str("");
        buffer << "tar -czvf " << directory << ".tar.gz "<< directory << "/ && rm -r " << directory << "/";
        silent_shell(buffer.str().c_str());
    }

    std::vector<std::string> sort_files_by_timestamp(std::string directory,std::string pattern)
    {
        std::set <fs::path> sort_by_name;
        for (fs::path p : fs::directory_iterator(directory))
        {
            if (p.extension() == pattern) 
            {
                sort_by_name.insert(p);
            }    
        }
        std::vector<std::string> file_list={};

        for (auto p : sort_by_name)
        {
            std::cout << p << std::endl;
            file_list.push_back(p);
        }
        
        return file_list;
    }
    bool file_exists_and_can_be_opened_and_is_not_empty(std::string filename)
    {
        if (! fs::exists(filename))
        {
            std::cout << "Input file not found: " << filename << std::endl;
            return false;
        }
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cout << "Input file cannot be opened for reading." << std::endl;
            return false;
        }
        file.close();
        if (! count_lines_in_file(filename))
        {
            std::cout << "Input file is empty: " << filename << std::endl;
            return false;
        }
        return true;
    }

    void ExitOnFailure(bool exit_if_false, std::string exit_message)
    {
        if (! exit_if_false)
        {
            std::cout << exit_message << std::endl;
            exit(0);
        }
    }

    int get_file_extension(std::string file_name)
    {
        std::string file_ext = fs::path(file_name).extension();
        file_ext = file_ext.substr(1,file_ext.size()-1);
        std::cout << "File extension is: " << file_ext << std::endl;

        std::vector <std::string> supported_filetypes = {"pdb","xyz","prmtop","mol2"};
        std::vector<std::string>::iterator ftype_pos = std::find(supported_filetypes.begin(),supported_filetypes.end(),file_ext);
        ExitOnFailure(ftype_pos != supported_filetypes.end(),"File extension is not supported.");
        int ftype_num =  ftype_pos - supported_filetypes.begin();
        return ftype_num;
    }
}
