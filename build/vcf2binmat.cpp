#include "utils.h"
#include "vcf2binmat.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>
#include <tuple>

/* Returns number of samples in vcf. */
unsigned int get_raw_count(const std::string vcf_path)
{
	std::ifstream vcf_file;
	vcf_file.open(vcf_path);
	if (!vcf_file.is_open())
		std::cerr << "Error opening file " << vcf_path << std::endl;
	
	const unsigned int offset = 9; // number of non-sample vcf columns
	std::string line;
	while (std::getline(vcf_file, line))
	{
		if (line.front() == '#') continue;
		auto split_line = splitString(line, '\t');
		return split_line.size() - offset;
	}
	return 0; // should never happen
}

/* Read vcf info file to select the indices needed to get a certain population or super population.
Call with empty string to leave it unspecified. */
std::vector<unsigned int> get_individual_indices(const std::string info_path,
	const std::string pop)
{
	std::cerr << "Retrieving relevant indices for vcf file ..." << std::endl;
	std::vector<unsigned int> result;
	std::ifstream info_file;
	info_file.open(info_path);
	if (!info_file.is_open())
		std::cerr << "Error opening file " << info_path << std::endl;
	
	std::string line, line_pop, line_super;
	unsigned int idx = 0;
	while (std::getline(info_file, line))
	{
		// header line
		if (line.find("gender") != std::string::npos)
			continue;
		
		idx += 1;
		
		if (pop.empty()) // no pop specified -> add all
			result.push_back(idx-1);
		else // only add on match
		{
			auto split_line = splitString(line, '\t');
			line_pop = split_line[1];
			line_super = split_line[2];
			
			if (line_pop.compare(pop) == 0
				|| line_super.compare(pop) == 0)
				result.push_back(idx-1);
		}
	}
	info_file.close();
	return result;
}

/* Returns the frequency of non-reference alleles or reference alleles
 * (i.e., non-0s and 0s) in given string, whichever is lower. */
double get_least_frequency(const std::string& s)
{
	auto zero_count = std::count(s.begin(), s.end(), '0');
	double freq = zero_count / (double) s.size();
	return std::min(freq, 1-freq);
}

/* Outputs genotypes of SNPs in VCF file to output file and records their
 * positions to an additional file. One line in the output file
 * corresponds to one variant site, one column to a haplotype.
 * Can be filtered by populations and minor allele frequency (MAF).
 * Returns: number of lines written, number of SNPs filter by MAF. */
std::tuple<unsigned int, unsigned int> vcf_to_binmat(
	const std::string vcf_path, const std::string out_path,
	const std::string info_path, const std::string pop, 
	const double maf_freq, const unsigned int maf_count)
{
	std::vector<unsigned int> indices;
	unsigned int individual_count = 0;
	if (!info_path.empty())
	{
		indices = get_individual_indices(info_path, pop);
		individual_count = indices.size();
		if (individual_count == 0)
		{
			std::cerr << "No individual found to match criteria!" << std::endl;
			return std::make_tuple(0, 0);
		}
	}
	// set individual count and indices to number of individuals in vcf if not set before
	if (individual_count == 0)
	{
		individual_count = get_raw_count(vcf_path);
		for (unsigned int i = 0; i < individual_count; i++)
			indices.push_back(i);
	}
	
	const double count_freq = maf_count / (double) (2*individual_count);
	double freq_threshold = std::max(maf_freq, count_freq);
	std::cerr << individual_count << " individuals  considered, MAF " << freq_threshold << std::endl;
	
	// open all the files
	std::ifstream vcf_file; // input vcf file
	std::ofstream out_file; // binary matrix output
	std::ofstream pos_file; // SNP positions output
	vcf_file.open(vcf_path);
	if (!vcf_file.is_open())
		std::cerr << "Error opening file " << vcf_path << std::endl;
	out_file.open(out_path, std::ofstream::out | std::ofstream::trunc);
	if (!out_file.is_open())
		std::cerr << "Error opening file " << out_path << std::endl;
	pos_file.open(out_path + ".positions", std::ofstream::out | std::ofstream::trunc);
	if (!pos_file.is_open())
		std::cerr << "Error opening file " << (out_path + ".positions") << std::endl;
	
	// the first nine entries of a vcf file line are not of interest
	unsigned int offset = 9, written_snps = 0, filtered_snps = 0;
	std::string line;
	while (std::getline(vcf_file, line))
	{
		// skip vcf header
		if (line.front() == '#')
			continue;
		auto split_line = splitString(line, '\t', individual_count+offset);
		
		// skip everything that is not a bi-allelic SNP
		if (split_line[3].size() != 1 || split_line[4].size() != 1)
			continue;
		
		std::string variant_line = "";
		variant_line.reserve(2*individual_count);
		// extract genotype info from valid line
		for (auto individual_id : indices)
		{
			// add offset to individual_id
			auto id = individual_id + offset;
			
			// example split_line[id]: "0|1"
			if (split_line[id].size() > 3)
				std::cerr << "Got overly long haplotype for individual "
					<< individual_id << ": " << split_line[id] << std::endl;

			try {
				variant_line += split_line[id].at(0);
				variant_line += split_line[id].at(2);
			}
			catch (const std::out_of_range& e)
			{
				std::cerr << e.what() << std::endl;
				std::cerr << "Thrown when trying to parse " << split_line[id] << "."
					<< " Possibly not diploid?" << std::endl;
				std::cerr << "SNP: " << split_line[2] << std::endl;
				std::cerr << "Individual: " << individual_id << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		
		// don't write if frequency is too low
		double freq = get_least_frequency(variant_line);
		if ((freq - freq_threshold) < (- 0.0001))
		// some leniency needed here due to rounding errors
		{
			filtered_snps++;
			continue;
		}
		
		// write to binary matrix and positions
		out_file << variant_line << "\n";
		pos_file << split_line[1] << "\n";
		written_snps++;
	}
	std::cerr << "\r" << "Wrote " << written_snps << " SNPs, filtered " << filtered_snps << std::endl;
	vcf_file.close();
	out_file.close();
	pos_file.close();
	return std::make_tuple(written_snps, filtered_snps);
}
