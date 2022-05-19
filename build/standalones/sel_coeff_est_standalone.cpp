#include "../utils.h"
#include "../sel_coeff_est.h"
#include <string>
#include <iostream>

class Config
{
public:
	Config(){};
	
	// output folder
	std::string out_file = "";
	
	// location of the vcf file
	std::string blocks_path = "";
	
	// location of the genetic map file (physical positions -> centimorgan)
	// e.g. http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/
	std::string genetic_map_path;
	// file with physical position for each SNP, ordered by position
	std::string pos_path = "";
	// file containing lookuptable for age thresholds to filter blocks
	std::string lookup_path = "";
	
	// effective population size
	unsigned int eff_pop_size = 10000;
	
	// maximal quantile for age filter
	double qmax = 0.01;
	// minimal selection coefficient we consider
	double smin = 0.005;
	
	// human chromosome number, will try to guess from file if not given
	unsigned int chr_number = 0;
	
	// if set, age and coalescence filter will be skipped
	bool skip_filters = false;
};

/* Parses cmd arguments into Config object. */
bool parse_command_args(Config &config, int argc, char* argv[])
{
	if (has_cmd_option(argv, argv+argc, "-h")
		|| !has_cmd_option(argv, argv+argc, "--out_file")
		|| !has_cmd_option(argv, argv+argc, "--blocks_path")
		|| !has_cmd_option(argv, argv+argc, "--pos_path")
		|| !has_cmd_option(argv, argv+argc, "--genetic_map_path")
		|| !has_cmd_option(argv, argv+argc, "--lookup_path"))
	{
		std::cerr << "Mandatory args: --out_file --blocks_path --pos_path --genetic_map_path --lookup_path" << std::endl;
		std::cerr << "Optional  args: --eff_pop_size --qmax --smin --chr_number --skip_filters" << std::endl;
		return false;
	}

	config.blocks_path = get_cmd_option(argv, argv+argc, "--blocks_path");
	config.out_file = get_cmd_option(argv, argv+argc, "--out_file");
	config.pos_path = get_cmd_option(argv, argv+argc, "--pos_path");
	config.genetic_map_path = get_cmd_option(argv, argv+argc, "--genetic_map_path");
	config.lookup_path = get_cmd_option(argv, argv+argc, "--lookup_path");
	
	if (has_cmd_option(argv, argv+argc, "--eff_pop_size"))
		config.eff_pop_size = std::stoul(get_cmd_option(argv, argv+argc, "--eff_pop_size"));
	if (has_cmd_option(argv, argv+argc, "--qmax"))
		config.qmax = std::stod(get_cmd_option(argv, argv+argc, "--qmax"));
	if (has_cmd_option(argv, argv+argc, "--smin"))
		config.smin = std::stod(get_cmd_option(argv, argv+argc, "--smin"));
	if (has_cmd_option(argv, argv+argc, "--chr_number"))
		config.chr_number = std::stoul(get_cmd_option(argv, argv+argc, "--chr_number"));
	config.skip_filters = has_cmd_option(argv, argv+argc, "--skip_filters");
	return true;
}

int main(int argc, char* argv[])
{
	Config config;
	if (!parse_command_args(config, argc, argv))
		return EXIT_FAILURE;
	
	estimate_selection_coeff(config.blocks_path, config.out_file,
		config.genetic_map_path, config.pos_path, config.lookup_path,
		config.eff_pop_size, config.skip_filters,
		config.qmax, config.smin, config.chr_number);
	return EXIT_SUCCESS;
}
