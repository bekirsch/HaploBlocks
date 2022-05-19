#include "utils.h"
#include "create_hist.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib> // std::system()
#include <sstream> // std::stringstream
#include <locale> // std::locale

// format big numbers (1000 -> 1.000)
template<class T>
std::string format_big_nr(T value)
{
	std::stringstream ss;
	ss.imbue(std::locale(""));
	ss << std::fixed << value;
	return ss.str();
}

// https://stackoverflow.com/questions/36288421/c-create-png-bitmap-from-array-of-numbers
// https://en.wikipedia.org/wiki/Netpbm_format
// use "convert image.pbm -scale 400x result.png" to convert
void print_pbm(const std::string in_path, const std::string out_path,
	const unsigned int resolution, unsigned int chr_number)
{
	std::ifstream stats_file;
	stats_file.open(in_path);
	if (!stats_file.is_open())
		std::cerr << "Error opening file " << in_path << std::endl;
	std::string line;
	std::vector<unsigned short> heights;
	
	// fill as much as wee need to if we have the chromosome nr
	if (chr_number == 0)
		chr_number = guess_chromosome_nr(out_path);
	if (chr_number > 0)
	{
		unsigned int needed_pixels = get_chromosome_size(chr_number) / resolution;
		while (needed_pixels+1000 > heights.size())
			heights.push_back(0);
	}
	
	// read stats file
	std::getline(stats_file, line); // skip first line
	while (std::getline(stats_file, line))
	{
		// read relevant info from line
		auto split_line = splitString(line, ',', 8);
		unsigned int bp_start = std::stoul(split_line[4]);
		unsigned int bp_end = std::stoul(split_line[5]);
		float s_hat = std::stof(split_line[8]);
		
		// convert info to pixel coordinates
		unsigned int pixel_start = bp_start / resolution;
		unsigned int pixel_end = bp_end / resolution;
		unsigned short pixel_height = s_hat * 5000;
		
		// update heights with current block
		while (pixel_end >= heights.size())
			heights.push_back(0);
		for (auto i = pixel_start; i <= pixel_end; i++)
			heights[i] = std::max(pixel_height, heights[i]);
	}
	stats_file.close();
	
	unsigned short max_height = *std::max_element(heights.begin(), heights.end());
	std::cerr << "Reading of input file done, writing image. "
		<< "Dimensions (width x height): " << heights.size() << " x " << max_height << std::endl;
	
	// shenanigans to make file name input robust ...
	std::string out_p(out_path);
	if (endswith(out_p, ".png"))
		out_p = out_p.substr(0, out_p.size()-4) + ".pbm";
	if (out_p.length() < 4)
		out_p += ".pbm";
	std::string out_png(out_p.substr(0, out_p.size()-4) + ".png");
	std::string tmp_png(out_p.substr(0, out_p.size()-4) + "_tmp.png");
	
	// clear output file
	std::ofstream out_file;
	out_file.open(out_p, std::ofstream::out | std::ofstream::trunc);
	out_file.close();
	// print histogram (top down)
	out_file.open(out_p, std::ofstream::out | std::ofstream::app);
	// header line
	out_file << "P1" << " " << heights.size() << " " << max_height << "\n";
	// actual content
	for (auto row = 0; row <= max_height; row++)
	{
		for (auto height : heights)
			out_file << (height >= (max_height-row));
		out_file << "\n";
	}
	out_file.close();
	
	// convert .pbm to .png
	std::string command;
	std::system(("convert " + out_p + " " + out_png).c_str());
	
	// create white image: "convert -size 800x800 xc:white white.png"
	command = "convert -size " + std::to_string(heights.size()) + "x30 xc:white " + tmp_png;
	std::system(command.c_str());
	
	// create annotation picture with physical positions
	// note that the annotated position is at the beginning of the annotation
	std::string plus2("+2 ");
	unsigned int annot_density = 200; // how many pixels are between numbers
	unsigned int numbers_needed = heights.size() / annot_density;
	for (size_t i = 0; i < numbers_needed; i++) 
	{
		std::string x_position = std::to_string(i*annot_density);
		std::string physical_position = format_big_nr(i*annot_density*resolution);
		command = ("convert " + tmp_png + " -gravity southwest -pointsize 16 -annotate +" + x_position + plus2 + physical_position + " " + tmp_png);
		std::system(command.c_str());
	}
	
	// combine images: "convert in-1.jpg in-2.jpg -append out.jpg"
	command = ("convert " + out_png + " " + tmp_png + " -append " + out_png);
	std::system(command.c_str());
	remove(out_p.c_str());
	remove(tmp_png.c_str());
}
