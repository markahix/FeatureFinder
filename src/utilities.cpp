
#include "utilities.h"

// LOGGING FUNCTIONS - Normal, Error, Debug-flagged.
void error_log(std::string message, int exit_code)
{
    std::cerr << message << std::endl;
    exit(exit_code);
}
void normal_log(std::string message)
{
    std::cout << message << std::endl;
}
void debug_log(std::string message)
{
    if (DEBUG)
    {
        std::cout << "DEBUG: " << message << std::endl;
    }
}
void splash_screen()
{
    std::cout << "###########################" << std::endl;
    std::cout << "# AGIMUS is in control... #" << std::endl;
    std::cout << "###########################" << std::endl;
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

// System Commands
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
void silent_shell(std::string cmd)
{
    silent_shell(cmd.c_str());
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
std::string GetSysResponse(std::string cmd)
{
    return GetSysResponse(cmd.c_str());
}
bool CheckProgAvailable(const char* program)
{
        std::string result;
        std::stringstream cmd;
        cmd.str("");
        cmd << "which " << program;
        debug_log(cmd.str());
        result=GetSysResponse(cmd.str().c_str());
        if (result.empty())
        {
            std::cout << "Missing program: " << program << std::endl;
            return false;
        }
        return true;
}
bool CheckProgAvailable(const char* program, const char* module)
{
        std::string result;
        std::stringstream cmd;
        cmd.str("");
        cmd << "module load " << module << "; which " << program;
        debug_log(cmd.str());
        result=GetSysResponse(cmd.str().c_str());
        if (result.empty())
        {
            std::cout << "Missing program: " << program << std::endl;
            return false;
        }
        return true;
}
bool CheckProgAvailable(std::string program)
{
    return CheckProgAvailable(program.c_str());
}
bool CheckProgAvailable(std::string program, std::string module)
{
    return CheckProgAvailable(program.c_str(),module.c_str());
}

// File Handling
bool validate_file(std::string filename)
{
    return validate_file(fs::absolute(filename));
}
bool validate_file(fs::path filename)
{
    if (!fs::exists(filename))
    {
        return false;
    }
    std::ifstream ifile(filename,std::ios::in);
    if (!ifile.is_open())
    {
        return false;
    }
    ifile.close();
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
    if (!fs::exists(inputfilename))
    {
        write_to_file(inputfilename,"");
    }
    std::ofstream outFile;
    outFile.open(inputfilename,std::ios::app);
    outFile << buffer;
    outFile.close();
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

// String Utilities
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
int is_empty(std::string s) //overload of above.
{
    return is_empty(s.c_str());
}
std::string string_between(std::string incoming, std::string first_delim, std::string second_delim)
{
    unsigned first = incoming.find(first_delim);
    if (first == std::string::npos)
    {
        return incoming;
    }
    incoming = incoming.substr(first+1,incoming.size() - first -1);

    unsigned last = incoming.find(second_delim);
    if (last == std::string::npos)
    {
        return incoming.substr(0, incoming.size());
    }
    return incoming.substr(0, last);
}
std::vector<std::string> split_string(std::string incoming, std::string delim)
{
    std::vector<std::string> chunks = {};
    unsigned int pos = incoming.find(delim);
    do
    {
        std::string piece = incoming.substr(0,pos);
        if (!is_empty(piece.c_str()))
        {
            chunks.push_back(piece);
        }       
        incoming = incoming.substr(pos+1,incoming.size()-pos-1);
        pos = incoming.find(delim);
    } while (pos != std::string::npos);
    chunks.push_back(incoming.substr(0,incoming.find_first_of(delim)));
    return chunks;
}



