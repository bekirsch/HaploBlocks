#include "../utils.h"
#include "../filter_blocks.h"
#include <string>
#include <iostream>

class Config
{
public:
	Config(){};
	
	// filtered csv output file
	std::string out_file = "";
	
	// bedGraph output file
	std::string bedgraph_file = "";
	
	// bed output file
	std::string bed_file = "";
	
	// whether bedGraph output should have overlaps or not
	bool overlap = false;
	
	// location of the input file
	std::string csv_path = "";
};

/* Parses cmd arguments into Config object. */
bool parse_command_args(Config &config, int argc, char* argv[])
{
	if (has_cmd_option(argv, argv+argc, "-h")
		|| !has_cmd_option(argv, argv+argc, "--out_file")
		|| !has_cmd_option(argv, argv+argc, "--csv_path"))
	{
		std::cerr << "Mandatory args: --out_file --csv_path" << std::endl;
		std::cerr << "Optional: --bedgraph_file if filtered csv should be turned into bedGraph." << std::endl;
		std::cerr << "Optional: --bed_file if filtered csv should be turned into bed." << std::endl;
		return false;
	}

	config.csv_path = get_cmd_option(argv, argv+argc, "--csv_path");
	config.out_file = get_cmd_option(argv, argv+argc, "--out_file");
	if (has_cmd_option(argv, argv+argc, "--bedgraph_file"))
		config.bedgraph_file = get_cmd_option(argv, argv+argc, "--bedgraph_file");
	if (has_cmd_option(argv, argv+argc, "--bed_file"))
		config.bed_file = get_cmd_option(argv, argv+argc, "--bed_file");
	config.overlap = has_cmd_option(argv, argv+argc, "--overlap");
	return true;
}

int main(int argc, char* argv[])
{
	Config config;
	if (!parse_command_args(config, argc, argv))
		return EXIT_FAILURE;
	
	filter_blocks(config.csv_path, config.out_file);
	if (config.bedgraph_file != "")
		csv2bedgraph(config.out_file, config.bedgraph_file, config.overlap);
	if (config.bed_file != "")
		csv2bed(config.out_file, config.bed_file, config.overlap);
	
	return EXIT_SUCCESS;
}
