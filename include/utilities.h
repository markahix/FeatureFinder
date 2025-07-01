#ifndef utilities_h
#define utilities_h

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <experimental/filesystem>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <locale>
#include <vector>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include <iomanip>
#include <ctime>
#include <set>

namespace fs = std::experimental::filesystem;

namespace utils
{
    void splash_screen();
    void silent_shell(const char* cmd);
    std::string GetSysResponse(const char* cmd);
    bool CheckProgAvailable(const char* program);
    void write_to_file(std::string inputfilename, std::string buffer);
    void append_to_file(std::string inputfilename, std::string buffer);
    bool CheckFileExists(std::string filename);
    std::string GetTimeAndDate();
    int is_empty(const char *s);
    std::string LastLineOfFile(std::string filename);
    int count_lines_in_file(std::string filename);
    std::string string_between(std::string incoming, std::string first_delim, std::string second_delim);
    void compress_and_delete(std::string directory);
    std::vector<std::string> sort_files_by_timestamp(std::string directory,std::string pattern);
    bool file_exists_and_can_be_opened_and_is_not_empty(std::string filename);
    void ExitOnFailure(bool exit_if_false, std::string exit_message);
    int get_file_extension(std::string file_name);
}


#endif