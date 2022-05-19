#include "../utils.h"
#include "../create_hist.h"
#include <string>
#include <iostream>

class Config
{
public:
	Config(){};
	
	// output file
	std::string out_file = "";
	
	// location of the input file
	std::string coeff_path = "";
	
	// resolution of the histogram output
	// a value of 1000 means that 1000 bases will correspond to a single pixel
	unsigned int resolution = 5000;
	
	// human chromosome number, will try to guess from file if not given
	unsigned int chr_number = 0;
};

/* Parses cmd arguments into Config object. */
bool parse_command_args(Config &config, int argc, char* argv[])
{
	if (has_cmd_option(argv, argv+argc, "-h")
		|| !has_cmd_option(argv, argv+argc, "--out_file")
		|| !has_cmd_option(argv, argv+argc, "--coeff_path"))
	{
		std::cerr << "Mandatory args: --out_file --coeff_path" << std::endl;
		std::cerr << "Optional  args: --resolution --chr_number" << std::endl;
		return false;
	}
	
	config.coeff_path = get_cmd_option(argv, argv+argc, "--coeff_path");
	config.out_file = get_cmd_option(argv, argv+argc, "--out_file");
	
	if (has_cmd_option(argv, argv+argc, "--resolution"))
		config.resolution = std::stoul(get_cmd_option(argv, argv+argc, "--resolution"));
	if (has_cmd_option(argv, argv+argc, "--chr_number"))
		config.chr_number = std::stoul(get_cmd_option(argv, argv+argc, "--chr_number"));
	return true;
}

int main(int argc, char* argv[])
{
	Config config;
	if (!parse_command_args(config, argc, argv))
		return EXIT_FAILURE;
	
	print_pbm(config.coeff_path, config.out_file, config.resolution);
	return EXIT_SUCCESS;
}
