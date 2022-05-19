#pragma once
#include <string>
#include <tuple>

std::tuple<unsigned int, unsigned int> vcf_to_binmat(
	const std::string vcf_path, const std::string out_path,
	const std::string info_path, const std::string pop,
	const double maf_freq = 0.0, const unsigned int maf_count = 0);
