#include "utils.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <sys/stat.h>

/* Splits input string at each occurence of char delim. */
std::vector<std::string> splitString(const std::string& input, const char delim, unsigned int expectedSize)
{
	std::string buffer{ "" };
	std::vector<std::string> result;
	result.reserve(expectedSize);
	for (auto c : input)
	{
		if (c != delim)
			buffer += c;
		else
		{
			result.push_back(buffer);
			buffer.clear();
		}
	}
	if (!buffer.empty()) // don't forget last part of string
		result.push_back(buffer);
	return result;
}

/* Check whether arguments contain given option. */
bool has_cmd_option(char **begin, char **end, std::string option)
{
	return std::find(begin, end, option) != end;
}

/* Returns argument for option if it exists. */
char* get_cmd_option(char **begin, char **end, std::string option)
{
	char **itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
		return *itr;
	else
		return 0;
}

/* Check whether given file exists. */
bool file_exists(const std::string& filename)
{
	struct stat buffer;   
	return (stat (filename.c_str(), &buffer) == 0); 
}

/* Check whether string a ends with string b */
bool endswith(const std::string a, const std::string b)
{
	// position where b would start if a ended with it
	int search_start = a.length() - b.length();
	if (search_start < 0) return false;
	return (a.compare(search_start, b.length(), b) == 0);
}

/* Tries to guess chromosome number from file name, returns 0 if not determined. */
unsigned int guess_chromosome_nr(const std::string& path)
{
	for (int i = 22; i > 0; i--) 
	{
		std::string to_search = "chr" + std::to_string(i);
		if (path.find(to_search) != std::string::npos)
			return i;
	}
	return 0;
}

// taken from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
// not in use here, as phase3 of the 1000 genomes project was mapped against GRCh37
std::map<unsigned int, unsigned int> chromosome_sizes_GRCh38 = {{1,248956422},{2,242193529},{3,198295559},{4,190214555},{5,181538259},
	{6,170805979},{7,159345973},{8,145138636},{9,138394717},{10,133797422},
	{11,135086622},{12,133275309},{13,114364328},{14,107043718},{15,101991189},
	{16,90338345},{17,83257441},{18,80373285},{19,58617616},{20,64444167},
	{21,46709983},{22,50818468}};

// taken from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
std::map<unsigned int, unsigned int> chromosome_sizes_GRCh37 = {{1,249250621},{2,243199373},{3,198022430},{4,191154276},{5,180915260},
	{6,171115067},{7,159138663},{8,146364022},{9,141213431},{10,135534747},
	{11,135006516},{12,133851895},{13,115169878},{14,107349540},{15,102531392},
	{16,90354753},{17,81195210},{18,78077248},{19,59128983},{20,63025520},
	{21,48129895},{22,51304566}};


/* Given a number 1 <= x <= 22, returns size of human chromosome x in bp. */
unsigned int get_chromosome_size(unsigned int number)
{
	if (number >= 1 && number <= 22)
		// use GRCh37 for 1000 genomes project phase 3!
		return chromosome_sizes_GRCh37[number];
	
	// invalid number
	std::cerr << "Invalid chromosome number: " << number << std::endl;
	return 0;
}
