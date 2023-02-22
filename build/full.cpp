#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <chrono> // stopping time
#include <iomanip>
#include <ctime>

#include "utils.h"
#include "vcf2binmat.h"
#include "pbwt.h"
#include "sel_coeff_est.h"
#include "filter_blocks.h"
#include "create_hist.h"

/* Stores all possible command line argument values. */
class Config
{
public:
	Config(){};
	
	// will remove tmp files after completion if true
	bool remove_tmp_files = false;
	
	// will try to use existing tmp files as input if false
	bool overwrite = false;
	
	// output folder
	std::string out_folder = "";
	
	// location of the vcf file
	std::string vcf_path = "";
	
	// location of the genetic map file (physical positions -> centimorgan)
	// e.g. http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/
	std::string genetic_map_path;
	// file containing lookuptable for age thresholds to filter blocks
	std::string lookup_path = "";
	
	// info file for individiuals in 1000 genomes project
	// integrated_call_samples_v3.20130502.ALL.panel
	std::string pop_info_path = "";
	// population / ancestry to select from vcf file
	std::string pop = "";
	
	// minimal allele count of reported blocks
	unsigned int min_allele_count = 0;
	// minimal sequence count of reported blocks
	unsigned int min_seq_count = 0;
	
	// effective population size
	unsigned int eff_pop_size = 10000;
	// maximal quantile for age filter
	double qmax = 0.01;
	// minimal selection coefficient we consider
	double smin = 0.005;
	// skip allele age and coalescence filter in sHat calculation if true
	bool skip_filters = false;
	
	// min minor allele frequency (MAF)
	double min_frequency = 0.0;
	// min absolute count of minor allele, converted to MAF once sample size is known
	unsigned int min_count = 0;
	
	// resolution of histogram, value of 1000 means that 1kbp corresponds to one pixel
	unsigned int resolution = 5000;
	
	// attributes for small report to log file at end (secret feature)
	bool do_report = false;
	std::string report();
};

/* Parses cmd arguments into Config object. */
bool parse_command_args(Config &config, int argc, char* argv[])
{
	if (has_cmd_option(argv, argv+argc, "-h")
		|| !has_cmd_option(argv, argv+argc, "--out_folder")
		|| !has_cmd_option(argv, argv+argc, "--vcf_path")
		|| !has_cmd_option(argv, argv+argc, "--genetic_map_path")
		|| !has_cmd_option(argv, argv+argc, "--lookup_path"))
	{
		std::cerr << "Mandatory args: --out_folder --vcf_path --genetic_map_path --lookup_path" << std::endl;
		std::cerr << "Optional: --remove --overwrite --pop_info_path --pop --min_allele_count --min_seq_count --MAF --MAF_count --eff_pop_size --qmax --smin --skip_filters --resolution" << std::endl;
		return false;
	}
	if (has_cmd_option(argv, argv+argc, "--pop")
		&& !has_cmd_option(argv, argv+argc, "--pop_info_path"))
	{
		std::cerr << "--pop_info_path is mandatory if --pop is given." << std::endl;
		return false;
	}
	
	config.vcf_path = get_cmd_option(argv, argv+argc, "--vcf_path");
	config.genetic_map_path = get_cmd_option(argv, argv+argc, "--genetic_map_path");
	config.lookup_path = get_cmd_option(argv, argv+argc, "--lookup_path");
	config.out_folder = get_cmd_option(argv, argv+argc, "--out_folder");
	// make sure out_folder actually is a folder
	if (strcmp(&config.out_folder.back(), "/") != 0)
		config.out_folder += "/";
	
	if (has_cmd_option(argv, argv+argc, "--remove"))
		config.remove_tmp_files = true;
	if (has_cmd_option(argv, argv+argc, "--overwrite"))
		config.overwrite = true;
	if (has_cmd_option(argv, argv+argc, "--pop_info_path"))
		config.pop_info_path = get_cmd_option(argv, argv+argc, "--pop_info_path");
	if (has_cmd_option(argv, argv+argc, "--pop"))
		config.pop = get_cmd_option(argv, argv+argc, "--pop");
	if (has_cmd_option(argv, argv+argc, "--min_allele_count"))
		config.min_allele_count = std::stoul(get_cmd_option(argv, argv+argc, "--min_allele_count"));
	if (has_cmd_option(argv, argv+argc, "--min_seq_count"))
		config.min_seq_count = std::stoul(get_cmd_option(argv, argv+argc, "--min_seq_count"));
	if (has_cmd_option(argv, argv+argc, "--MAF"))
		config.min_frequency = std::stod(get_cmd_option(argv, argv+argc, "--MAF"));
	if (has_cmd_option(argv, argv+argc, "--MAF_count"))
		config.min_count = std::stoul(get_cmd_option(argv, argv+argc, "--MAF_count"));
	if (has_cmd_option(argv, argv+argc, "--eff_pop_size"))
		config.eff_pop_size = std::stoul(get_cmd_option(argv, argv+argc, "--eff_pop_size"));
	if (has_cmd_option(argv, argv+argc, "--qmax"))
		config.qmax = std::stod(get_cmd_option(argv, argv+argc, "--qmax"));
	if (has_cmd_option(argv, argv+argc, "--smin"))
		config.smin = std::stod(get_cmd_option(argv, argv+argc, "--smin"));
	config.skip_filters = has_cmd_option(argv, argv+argc, "--skip_filters");
	if (has_cmd_option(argv, argv+argc, "--resolution"))
		config.resolution = std::stoul(get_cmd_option(argv, argv+argc, "--resolution"));
	config.do_report = has_cmd_option(argv, argv+argc, "--report");
	return true;
}

/* Prints time passed since last_time in seconds to stderr */
unsigned int passed_time(std::chrono::time_point<std::chrono::steady_clock> &last_time)
{
	auto current_time = std::chrono::steady_clock::now();
	auto diff_sec = std::chrono::duration_cast<std::chrono::seconds>(current_time - last_time).count();
	std::cerr << diff_sec << " seconds." << std::endl;
	return diff_sec;
}

void do_everything(Config &config)
{
	auto begin = std::chrono::steady_clock::now();
	auto last_step = std::chrono::steady_clock::now();
	
	
	// get base name of input file
	size_t found = config.vcf_path.find_last_of("/\\");
	std::string base_name = config.vcf_path.substr(found+1);
	
	// file IDs
	std::string file_id_short = "";
	if (!config.pop.empty())
		file_id_short += "_" + config.pop;
	if (config.min_frequency > 0)
		file_id_short += "_MAF_" + std::to_string(config.min_frequency);
	if (config.min_count > 0)
		file_id_short += "_MAC_" + std::to_string(config.min_count);
	
	// log single report line at the very end
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	std::ostringstream oss;
	oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
	std::string report = oss.str();
	report += ";" + base_name + file_id_short + ";";
	
	// call vcf2binmat
	std::string bm_file = config.out_folder + base_name + file_id_short + ".bm";
	std::string pos_file = bm_file + ".positions";
	if (config.overwrite || !file_exists(bm_file))
	{
		auto res = vcf_to_binmat(config.vcf_path, bm_file,
			config.pop_info_path, config.pop,
			config.min_frequency, config.min_count);
		std::cerr << "VCF parsing done, took ";
		auto diff = passed_time(last_step);
		report += "vcf:wrote:" + std::to_string(std::get<0>(res)) +
			",filtered:" + std::to_string(std::get<1>(res)) +
			",seconds:" + std::to_string(diff);
	}
	last_step = std::chrono::steady_clock::now();
	
	// call pbwt (haploblock calculation)
	std::string block_file = config.out_folder + base_name + file_id_short + ".blocks";
	if (config.overwrite || !file_exists(block_file))
	{
		auto res = calc_haploblocks(bm_file, block_file, config.min_allele_count, config.min_seq_count);
		std::cerr << "Calculation of haploblocks done, took ";
		auto diff = passed_time(last_step);
		report += ";pbwt:wrote:" + std::to_string(res) +
			",seconds:" + std::to_string(diff);
	}
	if (config.remove_tmp_files) remove(bm_file.c_str());
	last_step = std::chrono::steady_clock::now();
	
	// call sel_coeff_est
	std::string coeff_file = config.out_folder + base_name + file_id_short + ".sHat.csv";
	if (config.overwrite || !file_exists(coeff_file))
	{
		auto res = estimate_selection_coeff(block_file, coeff_file,
			config.genetic_map_path, pos_file, config.lookup_path,
			config.eff_pop_size, config.skip_filters,
			config.qmax, config.smin);
		std::cerr << "Calculation of selection coefficients done, took ";
		auto diff = passed_time(last_step);
		report += ";sel:wrote:" + std::to_string(std::get<0>(res)) +
			",filtered:" + std::to_string(std::get<1>(res)) +
			",seconds:" + std::to_string(diff);
	}
	if (config.remove_tmp_files) remove(block_file.c_str());
	if (config.remove_tmp_files) remove(pos_file.c_str());
	last_step = std::chrono::steady_clock::now();
	
	// call filter_blocks
	std::string coeff_file_filt = config.out_folder + base_name + file_id_short + "_filtered.sHat.csv";
	if (config.overwrite || !file_exists(coeff_file_filt))
	{
		auto res = filter_blocks(coeff_file, coeff_file_filt);
		std::cerr << "Filtering of selection coefficients done, took ";
		auto diff = passed_time(last_step);
		report += ";filter:wrote:" + std::to_string(res) +
			",seconds:" + std::to_string(diff);
	}
	if (config.remove_tmp_files) remove(coeff_file.c_str());
	last_step = std::chrono::steady_clock::now();
	
	// call csv2bed
	std::string bed_file = config.out_folder + base_name + file_id_short + ".bed";
	if (config.overwrite || !file_exists(bed_file))
	{
		csv2bed(coeff_file_filt, bed_file, false);
		std::cerr << "Filtered csv converted to bed, took ";
		passed_time(last_step);
	}
	last_step = std::chrono::steady_clock::now();
	
	// call csv2bedgraph
	std::string bedGraph_file = config.out_folder + base_name + file_id_short + ".bedGraph";
	if (config.overwrite || !file_exists(bedGraph_file))
	{
		csv2bedgraph(coeff_file_filt, bedGraph_file, false);
		std::cerr << "Filtered csv converted to bedgraph, took ";
		passed_time(last_step);
	}
	last_step = std::chrono::steady_clock::now();
	
	// call create_hist
	//std::string file_id_res = file_id_short + "_" + std::to_string(config.resolution);
	//std::string hist_file = config.out_folder + base_name + file_id_res + ".pbm";
	//std::string hist_file_png = config.out_folder + base_name + file_id_res + ".png";
	//if (config.overwrite || !file_exists(hist_file_png))
	//{
	//	print_pbm(coeff_file_filt, hist_file, config.resolution);
	//	std::cerr << "Histogram created, took ";
	//	passed_time(last_step);
	//}
	
	if (config.do_report)
	{
		std::ofstream outfile;
		outfile.open("full.log", std::ios_base::app);
		outfile << report << "\n";
		outfile.close();
	}
	
	std::cerr << "Program done, total time taken: ";
	passed_time(begin);
	return;
}

int main(int argc, char* argv[])
{
	Config config;
	if (!parse_command_args(config, argc, argv))
		return EXIT_FAILURE;
	
	do_everything(config);
	return EXIT_SUCCESS;
}
