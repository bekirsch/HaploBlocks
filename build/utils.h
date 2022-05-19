#pragma once
#include <vector>
#include <string>

/* Splits input string at each occurence of char delim. */
std::vector<std::string> splitString(const std::string& input, char delim, unsigned int expectedSize = 1);

/* Check whether arguments contain given option. */
bool has_cmd_option(char **begin, char **end, std::string option);

/* Returns argument for option if it exists. */
char* get_cmd_option(char **begin, char **end, std::string option);

/* Check whether given file exists. */
bool file_exists(const std::string& filename);

/* Check whether string a ends with string b */
bool endswith(const std::string a, const std::string b);

/* Tries to guess chromosome number from file name, returns 0 if not determined. */
unsigned int guess_chromosome_nr(const std::string& path);

/* Given a number 1 <= x <= 22, returns size of human chromosome x in bp */
unsigned int get_chromosome_size(unsigned int number);
