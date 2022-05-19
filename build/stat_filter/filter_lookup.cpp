#include "../utils.h"
#include "brent.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>	// std::max_element, std::max
#include <numeric>		// std::accumulate
#include <cfloat>		// DBL_EPSILON
#include <boost/math/distributions/gamma.hpp>


/** Config. **/
// effective population size
const unsigned int N_e_CONST = 20000;

/* Stores all possible command line argument values. */
class Config
{
public:
	Config(){};
	
	// quantile for generation age
	double q = 1.0e-02;
	
	// bounds and step for number of sequences k
	unsigned int min_k = 2;
	unsigned int max_k = 5008;
	
	// bounds and step for recombination fraction r
	double min_r = DBL_EPSILON;
	double max_r = 0.04;
	double step_r = 1e-9; // MINIMAL step, might be higher (see function next_r below)
};

/* Function for a smooth increase in k. */
unsigned int next_k(const unsigned int k, const unsigned int max_k)
{
	if (k >= max_k) return k + 1;
	unsigned int result = k * 1.05 + 1; // k * 1.05 rounds down, hence +1
	if (result > max_k) result = max_k; // make sure to not overshoot
	return result;
}

/* Function for a smooth increase in r. */
double next_r(const double r, const double step_r) 
{
	double constant = r + step_r;
	double smooth = r * 1.05;
	return std::max(constant, smooth);
}

/* Parses cmd arguments into Config object. */
bool parse_command_args(Config &config, int argc, char* argv[])
{
	if (has_cmd_option(argv, argv+argc, "-h"))
	{
		std::cerr << "Optional  args: --max_k --max_r --step_r --q" << std::endl;
		return false;
	}

	if (has_cmd_option(argv, argv+argc, "--max_k"))
		config.max_k = std::stoul(get_cmd_option(argv, argv+argc, "--max_k"));
	
	if (has_cmd_option(argv, argv+argc, "--max_r"))
		config.max_r = std::stof(get_cmd_option(argv, argv+argc, "--max_r"));
	if (has_cmd_option(argv, argv+argc, "--step_r"))
		config.step_r = std::stof(get_cmd_option(argv, argv+argc, "--step_r"));
	
	if (has_cmd_option(argv, argv+argc, "--q"))
		config.q = std::stof(get_cmd_option(argv, argv+argc, "--q"));
	
	// sanity checks to prevent infinite loops
	if (config.step_r <= DBL_EPSILON)
	{
		std::cerr << "Don't choose zero as step!" << std::endl;
		return false;
	}
	
	return true;
}

/** Statistical functions. **/
/* construct list of scale parameters for beta equation */
std::vector<double> beta_list(const unsigned int k, const double r)
{
	std::vector<double> result;
	result.reserve(k-1);
	for (size_t i = 2; i <= k; i++)
		result.push_back((2 * N_e_CONST) / (i * (i - 1 + 2 * N_e_CONST * r)));
	return result;
}

/* function to find root for as scale for gamma distribution */
double beta_equation(const double beta, const double mu, const std::vector<double> betas)
{
	double beta_sum = 0;
	for (auto b : betas)
		beta_sum += pow(b, 3) / pow(b + beta, 2);
	return mu / 2 - 2 * beta_sum;
}
/* Helper class for using brent.hpp  */
class betaEqFunctor : public brent::func_base
{
public:
	const double mu;
	const std::vector<double> betas;
	
	betaEqFunctor(const double &mu, const std::vector<double> &betas) :
		mu(mu), betas(betas) {}
	
	double operator() (double x)
	{
		return beta_equation(x, mu, betas);
	}
};

/* Approximate convolution via gamma approximation.
 * 
 * Returns beta_hat and writes mu to argument mu. */
double gammaprox(const unsigned int k, const double r, double &mu)
{
	std::vector<double> betas = beta_list(k, r);
	mu = std::accumulate(betas.begin(), betas.end(), 0.0f);
	double lower_bound = mu / betas.size(); // betas.size() = k - 1
	double upper_bound = *std::max_element(betas.begin(), betas.end());
	
	// search for argument that minimized the beta equation
	double beta_hat = 0;
	betaEqFunctor functor(mu, betas);
	double precision = sqrt(DBL_EPSILON);
	beta_hat = brent::zero(lower_bound, upper_bound, precision, functor);

	return beta_hat;
}

/* Given number of sequences in a block (k) and the recombination fraction
 * between the block's ends (r), calculate the quantile q for the age
 * of the block.
 * 
 * Using this as a filter, we shall only accept blocks younger than
 * the quantile returned by this function, as this suggests selection
 * as a driving factor rather than drift. */
double get_quantile(const unsigned int k, const double r, const double q)
{
	double mu = 0;
	double beta_hat = gammaprox(k, r, mu);
	double shape = mu / beta_hat;
	double scale = beta_hat;
	boost::math::gamma_distribution<double> distribution(shape, scale);
	return quantile(distribution, q);
}


/* Outputs lookup table for dimensions specified in config.
 * x axis has recombination fraction (r) increasing left to right
 * y axis has number of sequences (k) increasing top to bottom */
void write_lookup(const Config config)
{
	// calculate the upper bound for r to account for floating point errors
	double max_r = config.max_r + config.step_r / 2;
	
	// write header line
	for (auto r = config.min_r; r <= max_r; r = next_r(r, config.step_r))
		std::cout << "\t" << r;
	std::cout << std::endl;
	
	// write content
	for (auto k = config.min_k; k <= config.max_k; k = next_k(k, config.max_k))
	{
		std::cout << k;
		for (auto r = config.min_r; r <= max_r; r = next_r(r, config.step_r))
			std::cout << "\t" << get_quantile(k, r, config.q);
		std::cout << std::endl;
	}
}

int main(int argc, char *argv[])
{
	Config config;
	if (!parse_command_args(config, argc, argv))
		return EXIT_FAILURE;
	
	write_lookup(config);
	return EXIT_SUCCESS;
}
