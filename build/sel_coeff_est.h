#pragma once
#include <string>
#include <tuple>

std::tuple<unsigned int, unsigned int> estimate_selection_coeff(
	const std::string in_path, const std::string out_path,
	const std::string recomb_map_path, const std::string pos_list_path,
	const std::string lookup_path, const unsigned int eff_pop_size,
	const bool skip_filters, const double qmax, const double smin,
	unsigned int chr_number = 0);
