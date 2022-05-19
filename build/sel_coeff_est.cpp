#include "utils.h"
#include "sel_coeff_est.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <utility>
#include <cmath>
#include <climits>
#include <algorithm>
#include <map>
#include <numeric>      // std::adjacent_difference, std::accumulate

// example input: "1013958 (23549-24439,3014:1138)"
// "weight (start-end,witness:seq_count)"
void extract_block_info(const std::string &input,
	unsigned int &weight, unsigned int &start, unsigned int &end, unsigned int &witness, unsigned int &seq_count)
{
	std::string buffer{""};
	unsigned int idx = 0;
	while (input[idx] != ' ')
	{
		buffer += input[idx];
		idx++;
	}
	weight = std::stoul(buffer);
	buffer.clear();
	idx += 2; // skip open parenthesis
	while (input[idx] != '-')
	{
		buffer += input[idx];
		idx++;
	}
	start = std::stoul(buffer);
	buffer.clear();
	idx++;
	while (input[idx] != ',')
	{
		buffer += input[idx];
		idx++;
	}
	end = std::stoul(buffer);
	buffer.clear();
	idx++;
	while (input[idx] != ':')
	{
		buffer += input[idx];
		idx++;
	}
	witness = std::stoul(buffer);
	buffer.clear();
	idx++;
	while (input[idx] != ')')
	{
		buffer += input[idx];
		idx++;
	}
	seq_count = std::stoul(buffer);
}

/** Statistical functions. **/
/* Haldane map - centimorgan to recombination fraction*/
inline double cM_to_recomb(const double cM)
{
	return (1 - std::exp((-2 * cM) / 100)) / 2;
}

/* Estimate age of block given selection coefficient s, number of 
 * sequences in block n, total population size N,
 * and effective population size N_e. */
inline double t_given_s_n_N(const double s, const double n, const double N, const double N_e)
{
	// reminder: n/N = p, the frequency of the block in the sample
	return (1 / s) * std::log(((n/N) * (1 - 1/N_e)) / ((1/N_e) * (1 - n/N)));
}

/* Age estimation for frequency filtering. 
 * y is frequency, n sample size, q quantile. */
inline double t_freq(const double y, const double n, const double q)
{
	return (2 * std::log(1-y)) / (std::log(q) + std::log(1-y)) - 2/n;
}

double adaptive_quantile(const double n, const double N, const double s_min, const double N_e, const double q_max=0.01)
{
	double p = n / N; // block frequency
	double t_q = t_given_s_n_N(s_min, n, N, 2*N_e);
	t_q /= 2*N_e; // (we need t in units of 2*N_e generations)
	double qexp = N / (1 + N * t_q / 2) -1;
	double qs = std::pow(1-p, qexp);
	return std::min(q_max, std::max(0.0001, qs));
}

/* Calculate selection coefficient on a block containing block_size
 * sequences and having bounds with genetic distance delta_cM. */
double sHat(const int nr_chroms, const int block_size, const double r, const double delta_r, const unsigned int N_e)
{
	if (nr_chroms < block_size)
	{
		std::cerr << "ERROR: Cannot have more sequences in a block than there are sequences total." << std::endl;
		return 0;
	}
	if (r < 0)
	{
		std::cerr << "ERROR: recombination fraction cannot be negative." << std::endl;
		return 0;
	}
	
	const double y_0 = 1.0 / (2*N_e);
	const double y_t = (double) block_size / nr_chroms;
	
	const double factor1 = delta_r / std::log(r / (delta_r + r));
	const double factor2 = std::log(y_0 / y_t);
	return factor1 * factor2;
}

/** Filtering of blocks (selection inference). **/
struct LookUpTable
{
	std::vector<double> r_idx;
	std::vector<unsigned int> k_idx;
	std::vector<std::vector<double>> t_values;
	
	LookUpTable(std::vector<double> r_idx, std::vector<unsigned int> k_idx, std::vector<std::vector<double>> t_values) 
		: r_idx(r_idx), k_idx(k_idx), t_values(t_values) {};
	
	double lookup(const double r, const unsigned int k) const;
};

/* Look up table value in table, interpolating if necessary. */
double LookUpTable::lookup(const double r, const unsigned int k) const
{
	unsigned int r1, r2, k1, k2; // indices
	double r1v, r2v, k1v, k2v; //values
	
	r2 = std::distance(r_idx.begin(), std::lower_bound(r_idx.begin(), r_idx.end(), r));
	k2 = std::distance(k_idx.begin(), std::lower_bound(k_idx.begin(), k_idx.end(), k));
	r1 = (r2 == 0) ? 0 : (r2 - 1);
	k1 = (k2 == 0) ? 0 : (k2 - 1);
	// checks for lower OOB
	if (r1 == r2) r2++;
	if (k1 == k2) k2++;
	
	r1v = r_idx[r1]; r2v = r_idx[r2];
	k1v = k_idx[k1]; k2v = k_idx[k2];
	
	double q11, q12, q21, q22;
	q11 = t_values[k1][r1];
	q12 = t_values[k1][r2];
	q21 = t_values[k2][r1];
	q22 = t_values[k2][r2];
	
	// see https://en.wikipedia.org/wiki/Bilinear_interpolation
	double result = (q11 * (k2v-k) * (r2v-r) 
					+q21 * (k-k1v) * (r2v-r)
					+q12 * (k2v-k) * (r-r1v)
					+q22 * (k-k1v) * (r-r1v))
					/ ((k2v-k1v) * (r2v-r1v));
	
	return result;
}

/* Reads lookuptable from file and stores it as vector of vectors.
 * recombination fraction increases with columns.
 * number of sequences increases with rows. */
LookUpTable read_lookuptable(const std::string lookuptable_path)
{
	// open inputfile
	std::ifstream lookup;
	lookup.open(lookuptable_path);
	if (!lookup.is_open())
		std::cerr << "Error opening file " << lookuptable_path << std::endl;
	std::string line;
	
	std::vector<double> r_idx;
	std::vector<unsigned int> k_idx;
	std::vector<std::vector<double>> t_values;
	
	// first line, contains values for recombination fraction
	std::getline(lookup, line);
	for (auto val : splitString(line, '\t'))
	{
		if (val.size() < 1) continue;
		r_idx.push_back(std::stod(val));
	}
	// add extra column for out of bounds
	r_idx.push_back(1);
	auto nr_columns = r_idx.size();
	
	// remaining lines, filling k_idx and t_values
	while(std::getline(lookup, line))
	{
		if (line.size() < 1) continue;
		auto split_line = splitString(line, '\t', nr_columns+1);
		// first entry is k index
		k_idx.push_back(std::stoul(split_line[0]));
		
		// remaining entries are t values
		std::vector<double> t_row;
		t_row.reserve(nr_columns);
		for (size_t i = 1; i < split_line.size(); i++)
			t_row.push_back(std::stod(split_line[i]));
		// extra column for out of bound
		t_row.push_back(t_row.back());
		t_values.push_back(t_row);
	}
	
	// extra line for out of bound
	k_idx.push_back(std::numeric_limits<int>::max());
	t_values.push_back(t_values.back());
	
	return LookUpTable(r_idx, k_idx, t_values);
}

// given ordered xs, returns index at which x would need to be inserted to maintain order
// since idx only increases for ordered files, we remember the last one for efficiency
unsigned int last_idx = 0;
unsigned int get_index(const unsigned int x, const std::vector<unsigned int> &xs)
{
	unsigned int idx = last_idx;
	while (idx < xs.size())
	{
		if (x < xs[idx]) break;
		idx++;
	}
	last_idx = idx;
	return idx;
}

double linear_interpolation(const unsigned int x, const std::vector<unsigned int> &xs, const std::vector<double> &ys)
{
	auto idx = get_index(x, xs);
	auto x_0 = xs[idx - 1], x_1 = xs[idx];
	auto y_0 = ys[idx - 1], y_1 = ys[idx];
	double y = y_0 + (x - x_0) * (y_1 - y_0) / (x_1 - x_0);
	return y;
}

void parse_recombination_map(const std::string in_path, 
	std::vector<unsigned int> &base_positions, std::vector<double> &centimorgan)
{
	// first entries are left border
	centimorgan.push_back(0.0);
	base_positions.push_back(0);

	std::ifstream in_file;
	in_file.open(in_path);
	if (!in_file.is_open())
		std::cerr << "Error opening file " << in_path << std::endl;
	std::string line;
	while (std::getline(in_file, line))
	{
		if (line.front() == '#') continue;
		auto split_line = splitString(line, '\t', 4);
		if (split_line.size() != 4) continue;
		centimorgan.push_back(std::stod(split_line[2]));
		base_positions.push_back(std::stoul(split_line[3]));
	}
	in_file.close();
	// arbitrary right border to avoid extrapolation errors
	centimorgan.push_back(centimorgan.back());
	base_positions.push_back(UINT_MAX);
}

// builds map from vcf index to base position and centimorgen position
// e.g. index_to_cm[0] will hold the pair (x,y) where x is the base position of the first SNP and y its cM
void parse_physical_positions(const std::string in_path, std::vector<std::pair<unsigned int, double>> &index_to_cm,
	const std::vector<unsigned int> &base_positions, const std::vector<double> &centimorgan)
{
	std::ifstream in_file;
	in_file.open(in_path);
	if (!in_file.is_open())
		std::cerr << "Error opening file " << in_path << std::endl;
	std::string line;
	while (std::getline(in_file, line))
	{
		if (line.front() == '#') continue;
		if (line.size() < 1) continue;
		unsigned int position = std::stoul(line);
		auto cM = linear_interpolation(position, base_positions, centimorgan);
		index_to_cm.push_back(std::make_pair(position, cM));
	}
	in_file.close();
}

// calculates delta_r as the mean recombination fraction between all pairs of contigious SNPs 
double calc_deltar(const std::vector<std::pair<unsigned int, double>> &index_to_cm)
{
	// extract centimorgan values only
	std::vector<double> cm_continuous_SNPs;
	cm_continuous_SNPs.reserve(index_to_cm.size());
	for (auto pair : index_to_cm)
		cm_continuous_SNPs.push_back(pair.second);
	
	// write distance to next position (= pairwise cm distances) to cm_continuous_SNPs
	std::adjacent_difference(cm_continuous_SNPs.begin(), cm_continuous_SNPs.end(), cm_continuous_SNPs.begin());
	
	// calculate r_frac on each pairwise distance
	std::transform(cm_continuous_SNPs.begin(), cm_continuous_SNPs.end(), cm_continuous_SNPs.begin(), cM_to_recomb);
	
	// +-1 in the next to lines because the first value of adj_diff should be ignored
	double sum = std::accumulate(cm_continuous_SNPs.begin()+1, cm_continuous_SNPs.end(), 0.0);
	
	// return average
	return sum / (cm_continuous_SNPs.size()-1);
}

// it's super important that map_path and pos_path files are ordered by start position!
std::tuple<unsigned int, unsigned int> estimate_selection_coeff(
	const std::string in_path, const std::string out_path,
	const std::string recomb_map_path, const std::string pos_path,
	const std::string lookup_path, const unsigned int eff_pop_size,
	const bool skip_filters, const double qmax, const double smin,
	unsigned int chr_number)
{
	std::vector<double> cMs; std::vector<unsigned int> bPs;
	parse_recombination_map(recomb_map_path, bPs, cMs);
	std::vector<std::pair<unsigned int, double>> index_to_cm;
	parse_physical_positions(pos_path, index_to_cm, bPs, cMs);
	LookUpTable lookuptable = read_lookuptable(lookup_path);
	
	unsigned int nr_chroms = 5008; // assumption if not given in file, else set below
	unsigned int nr_alleles = 42; // assumption if not given in file, else set below
	std::ifstream in_file;
	in_file.open(in_path);
	if (!in_file.is_open())
		std::cerr << "Error opening file " << in_path << std::endl;
	std::string line;
	unsigned int weight, start, end, witness, count;
	unsigned int bp_start, bp_end;
	double cm_start, cm_end;
	
	if (chr_number == 0)
		chr_number = guess_chromosome_nr(in_path);
	const double delta_r = calc_deltar(index_to_cm);
	
	// clear & open output file
	std::ofstream out_file;
	out_file.open(out_path, std::ofstream::out | std::ofstream::trunc);
	unsigned int output_count = 0, age_filter_count = 0, coa_filter_count = 0;
	out_file << "witness,start,end,size,bp start,bp end,r frac,t [gen],sHat,freq\n";
	while (std::getline(in_file, line))
	{
		if (line.size() < 2) continue;
		
		if (line.size() > 0 && line.front() == '#')
		{ // header line containing original binary matrix dimension
			line.erase(0,1);
			auto infos = splitString(line, ',', 2);
			nr_chroms = std::stoul(infos[0]);
			if (infos.size() > 1)
				nr_alleles = std::stoul(infos[1]);
			continue;
		}
		
		extract_block_info(line, weight, start, end, witness, count);
		
		try
		{
			std::tie(bp_start, cm_start) = index_to_cm[start];
			std::tie(bp_end, cm_end) = index_to_cm[end];
		}
		catch (const std::exception& e)
		{
			std::cerr << "Exception when trying to access index_to_cm: " << e.what() << std::endl;
			std::cerr << "start=" << start << ", end=" << end << ", array size=" << index_to_cm.size() << std::endl;
			std::tie(bp_end, cm_end) = index_to_cm.back();
		}
		double delta_cM = cm_end - cm_start;
		if (delta_cM == 0) continue; // will have sHat == 0
		const double r = cM_to_recomb(delta_cM);
		double s_hat = sHat(nr_chroms, count, r, delta_r, eff_pop_size);
		
		double t_2 = t_given_s_n_N(s_hat, count, nr_chroms, eff_pop_size);
		if (!skip_filters)
		{ // coalescence filter
			double t_threshold = lookuptable.lookup(r, count);
			if (t_2 > t_threshold) 
			{
				coa_filter_count++;
				continue;
			}
		}
		
		double y = count / (double) nr_chroms;
		if (!skip_filters)
		{ // age filter
			double t_1 = t_given_s_n_N(s_hat, count, nr_chroms, 2*eff_pop_size);
			double q = adaptive_quantile(count, nr_chroms, smin, eff_pop_size, qmax);
			double t_threshold = t_freq(y, nr_alleles, q);
			if (t_1 > 2*eff_pop_size*t_threshold)
			{
				age_filter_count++;
				continue;
			}
		}
		
		// MAF in percent, rounded to two decimal points
		double maf = round(std::min(y, 1-y) * 10000.0) / 100.0;
		// output block
		out_file << witness << "," << start << "," << end  << "," <<
			count << "," << bp_start << "," << bp_end << "," << r << "," <<
			t_2 << "," << s_hat << "," << maf << "\n";
		output_count++;
	}
	unsigned int total_filter_count = age_filter_count + coa_filter_count;
	std::cerr << "\r" << output_count << " selection coefficient estimates written, "
		<< total_filter_count << " blocks filtered. (" << coa_filter_count
		<< " by coalescence, " << age_filter_count << " remaining by age)" << std::endl;
	out_file.close();
	in_file.close();
	return std::make_tuple(output_count, total_filter_count);
}
