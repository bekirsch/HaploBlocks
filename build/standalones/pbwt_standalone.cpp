#include "../utils.h"
#include "../pbwt.h"
#include <string>
#include <iostream>

class Config
{
public:
	Config(){};
	
	// output folder
	std::string out_file = "";
	
	// location of the vcf file
	std::string block_path = "";
	
	// minimal allele count of reported blocks
	unsigned int min_allele_count = 0;
	// minimal sequence count of reported blocks
	unsigned int min_seq_count = 0;
};

bool parse_command_args(Config &config, int argc, char* argv[])
{
	if (has_cmd_option(argv, argv+argc, "-h")
		|| !has_cmd_option(argv, argv+argc, "--out_file")
		|| !has_cmd_option(argv, argv+argc, "--block_path"))
	{
		std::cerr << "Mandatory args: --out_file --block_path" << std::endl;
		std::cerr << "Optional args: --min_allele_count --min_seq_count" << std::endl;
		return false;
	}
	
	config.block_path = get_cmd_option(argv, argv+argc, "--block_path");
	config.out_file = get_cmd_option(argv, argv+argc, "--out_file");
	
	if (has_cmd_option(argv, argv+argc, "--min_allele_count"))
		config.min_allele_count = std::stoul(get_cmd_option(argv, argv+argc, "--min_allele_count"));
	if (has_cmd_option(argv, argv+argc, "--min_seq_count"))
		config.min_seq_count = std::stoul(get_cmd_option(argv, argv+argc, "--min_seq_count"));
	return true;
}

int main(int argc, char* argv[])
{
	Config config;
	if (!parse_command_args(config, argc, argv))
		return EXIT_FAILURE;
	
	calc_haploblocks(config.block_path, config.out_file, config.min_allele_count, config.min_seq_count);
	return EXIT_SUCCESS;
}
