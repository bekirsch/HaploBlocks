#include "../utils.h"
#include "../vcf2binmat.h"
#include <string>
#include <iostream>

/* Stores all possible command line argument values.
 * Optional arguments have standards.
 * */
class Config
{
public:
	Config(){};
	
	// location of the vcf file
	std::string vcf_path;
	
	// output file
	std::string out_file;
	
	// file with physical position for each SNP, ordered by position
	std::string pos_path;
	
	// minor allele frequency (MAF)
	double min_frequency = 0.0;
	// min absolute count of minor allele
	// converted to MAF in code once population size is known
	unsigned int min_count = 0;
	
	// info file for individiuals in 1000 genomes project (e.g. integrated_call_samples_v3.20130502.ALL.panel)
	std::string pop_info_path;
	// population / super_population to select from vcf file
	std::string pop = "";
};

bool parse_command_args(Config &config, int argc, char* argv[])
{
	if (has_cmd_option(argv, argv+argc, "-h")
		|| !has_cmd_option(argv, argv+argc, "--out_file")
		|| !has_cmd_option(argv, argv+argc, "--vcf_path"))
	{
		std::cerr << "Mandatory args: --out_file --vcf_path" << std::endl;
		std::cerr << "Optional  args: --pop_info_path --pop --MAF --MAF_count" << std::endl;
		return false;
	}
	if (has_cmd_option(argv, argv+argc, "--pop")
		&& !has_cmd_option(argv, argv+argc, "--pop_info_path"))
	{
		std::cerr << "--pop_info_path is mandatory if --pop is given." << std::endl;
		return false;
	}

	config.vcf_path = get_cmd_option(argv, argv+argc, "--vcf_path");
	config.out_file = get_cmd_option(argv, argv+argc, "--out_file");
	
	if (has_cmd_option(argv, argv+argc, "--pop_info_path"))
		config.pop_info_path = get_cmd_option(argv, argv+argc, "--pop_info_path");
	if (has_cmd_option(argv, argv+argc, "--pop"))
		config.pop = get_cmd_option(argv, argv+argc, "--pop");
	if (has_cmd_option(argv, argv+argc, "--MAF"))
		config.min_frequency = std::stod(get_cmd_option(argv, argv+argc, "--MAF"));
	if (has_cmd_option(argv, argv+argc, "--MAF_count"))
		config.min_count = std::stoul(get_cmd_option(argv, argv+argc, "--MAF_count"));
		
	return true;
}

int main(int argc, char* argv[])
{
	Config config;
	if (!parse_command_args(config, argc, argv))
		return EXIT_FAILURE;
	
	vcf_to_binmat(config.vcf_path, config.out_file,
		config.pop_info_path, config.pop,
		config.min_frequency, config.min_count);
	return EXIT_SUCCESS;
}
