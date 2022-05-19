#include "utils.h"
#include "filter_blocks.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include <memory>
#include <algorithm>

/* Struct containing info related to a haploblock. */
struct Block
{
	unsigned int witness;
	unsigned int snp_start;
	unsigned int snp_end;
	unsigned int size;
	unsigned int bp_start;
	unsigned int bp_end;
	double r_frac;
	double t_gen;
	double sHat;
	double frequency;
	
	Block(unsigned int bp_start, unsigned int bp_end, float sHat) :
		bp_start(bp_start), bp_end(bp_end), sHat(sHat)
		{}
	
	Block(std::string line);
	
	std::string to_csv_line();
	std::string to_bed_line(unsigned int chromNumber,
		unsigned int chromStart, unsigned int chromEnd);
	std::string bedoutput(unsigned int chromNumber,
		unsigned int chromStart, unsigned int chromEnd, bool bedgraph);
	
	
	friend std::ostream& operator<<(std::ostream& os, const Block& b)
	{
		os << b.bp_start << ' ' << b.bp_end << ' ' << b.sHat;
		return os;
	}
};

/* Construct Block from line in .csv. */
Block::Block(std::string line)
{
	auto split_line = splitString(line, ',', 10);
	witness = std::stoul(split_line[0]);
	snp_start = std::stoul(split_line[1]);
	snp_end = std::stoul(split_line[2]);
	size = std::stoul(split_line[3]);
	bp_start = std::stoul(split_line[4]);
	bp_end = std::stoul(split_line[5]);
	r_frac = std::stod(split_line[6]);
	t_gen = std::stod(split_line[7]);
	sHat = std::stod(split_line[8]);
	frequency = std::stod(split_line[9]);
}

/* Turns block into csv line. */
std::string Block::to_csv_line()
{
	std::ostringstream oss;
	oss << witness << "," << snp_start << "," << snp_end << "," << size 
		<< "," << bp_start << "," << bp_end << "," << r_frac << ","
		<< t_gen << "," << sHat << "," << frequency;
	return oss.str();
}

/* Turns block into bed line except for the leading chr field.
 * see autoSQL "haploblocks_bigbed.as". */
std::string Block::to_bed_line(unsigned int chromNumber, unsigned int chromStart, unsigned int chromEnd)
{
	std::ostringstream oss;
	oss << "chr" << chromNumber << " " <<
		chromStart << " " << chromEnd << " " << frequency << " " <<
		size << " " << bp_start << " " << bp_end << " " <<
		(snp_end - snp_start + 1) << " " << witness << " " <<
		r_frac << " " << sHat << " " << t_gen;
	return oss.str();
}

std::string Block::bedoutput(unsigned int chromNumber, unsigned int chromStart, unsigned int chromEnd, bool bedgraph)
{
	if (bedgraph)
	{
		std::ostringstream oss;
		oss << "chr" << chromNumber << " " <<
			chromStart << " " << chromEnd << " " << sHat;
		return oss.str();
	}
	else
		return to_bed_line(chromNumber, chromStart, chromEnd);
}

/* Filters Histogram blocks from input file such that the output
 * does not contain any block that is completely covered by other
 * blocks, as they are redundant for the output.
 * Also removes block with height (sHat) 0.
 * 
 * This implies that a input file and the filtered output should
 * produce the same Histogram and also a visually identical bedGraph
 * for windowingFunction=maximum.
 * 
 * Returns number of blocks left after filtering.
 * */
unsigned int filter_blocks(const std::string in_path, const std::string out_path)
{
	std::ifstream stats_file;
	stats_file.open(in_path);
	if (!stats_file.is_open())
	{
		std::cerr << "Error opening file " << in_path << std::endl;
		return 0;
	}
	
	// *max_blocks[i] = block with highest shat at SNP i (or nullptr if no block at SNP)
	std::vector<std::shared_ptr<Block>> max_blocks;
	unsigned long count = 0;
	
	// read input file, only store max blocks
	std::string line;
	std::getline(stats_file, line); // skip first line
	while (std::getline(stats_file, line))
	{
		count++;
		std::shared_ptr<Block> ptr(new Block(line));
		if (ptr->sHat == 0) continue;
		
		for (auto i = ptr->snp_start; i <= ptr->snp_end; i++)
		{
			// SNP position not seen yet, cannot be compared
			if (i >= max_blocks.size()) break;
			
			float current_max = max_blocks[i] == nullptr ? 0.0 : max_blocks[i]->sHat;
			if (ptr->sHat > current_max) // new max value -> save
				max_blocks[i] = ptr;
		}
		
		// snp_start is beyond so far observed max SNP position (gap between)
		while (max_blocks.size() < ptr->snp_start)
			max_blocks.push_back(nullptr);
		
		// if vectors need to be extended, there were no blocks at those SNPs before
		while (max_blocks.size() <= ptr->snp_end)
			max_blocks.push_back(ptr);
	}
	stats_file.close();
	
	std::set<std::shared_ptr<Block>> unique_blocks(max_blocks.begin(), max_blocks.end());
	unique_blocks.erase(nullptr); // remove nullptr (not actually a Block)
	std::cout << "Input had " << count << " blocks, " << unique_blocks.size() << " left after filtering." << std::endl;
	
	// write unique blocks
	if (out_path.size() > 0)
	{
		std::ofstream out_file;
		out_file.open(out_path, std::ofstream::out | std::ofstream::trunc);
		if (!out_file.is_open())
		{
			std::cerr << "Error opening file " << out_path << std::endl;
			return unique_blocks.size();
		}
		out_file << "witness,start,end,size,bp start,bp end,cm delta,t [gen],sHat,freq\n";
		for (const auto& hp : unique_blocks)
			out_file << hp->to_csv_line() << '\n';
		out_file.close();
	}
	return unique_blocks.size();
}


/* Truncates input blocks such that there is no overlap and outputs them
 * sorted by their starting position.
 * Note that only base pair positions are affected.
 * If graphout prints bedGraph lines, otherwise bed lines. */
void print_without_overlap(std::vector<Block>& blocks, std::ofstream& out_file,
	const unsigned int chr_nr, bool graphout=false)
{
	if (blocks.size() < 2) return; // run with overlap = true ...
	// init here is similar to filter_blocks
	std::vector<Block*> max_blocks;
	
	for (size_t i = 0; i < blocks.size(); i++)
	{
		for (auto snp_pos = blocks[i].snp_start; snp_pos <= blocks[i].snp_end; snp_pos++)
		{
			// SNP position not seen yet, cannot be compared
			if (snp_pos >= max_blocks.size()) break;
			
			float current_max = max_blocks[snp_pos] == nullptr ? 0.0 : max_blocks[snp_pos]->sHat;
			if (blocks[i].sHat > current_max) // new max value -> save
				max_blocks[snp_pos] = &blocks[i];
		}
		
		// snp_start is beyond so far observed max SNP position (gap between)
		while (max_blocks.size() < blocks[i].snp_start)
			max_blocks.push_back(nullptr);
		
		// if vectors need to be extended, there were no blocks at those SNPs before
		while (max_blocks.size() <= blocks[i].snp_end)
			max_blocks.push_back(&(blocks[i]));
	}
	
	Block* prev_ptr = nullptr;
	unsigned int block_start = 0;
	for (auto ptr : max_blocks)
	{
		// nothing to do if both are the same
		if (prev_ptr == ptr) continue;
		
		// if prev_ptr was nullptr, we do not report it
		if (prev_ptr == nullptr)
		{
			prev_ptr = ptr;
			block_start = ptr->bp_start;
			continue;
		}
		
		// if new block is nullptr, simply report the previous block
		if (ptr == nullptr)
			out_file << prev_ptr->bedoutput(chr_nr, block_start, prev_ptr->bp_end, graphout) << "\n";
		else
		// two neighboring blocks, edit start / end to avoid overlap
		{
			// neighboring blocks with equal sHat may be redundant
			if (prev_ptr->sHat == ptr->sHat && prev_ptr->bp_end >= ptr->bp_end)
				continue;
			
			bool new_block_lower = (ptr->sHat < prev_ptr->sHat);
			unsigned int block_end = 0;
			if (new_block_lower)
				block_end = prev_ptr->bp_end;
			else
				block_end = ptr->bp_start - 1;
			
			/* There can be cases where neighboring blocks have equal sHat.
			 * These corner cases cause some weird behavior, which is fixed
			 * by the block_start < block_end condition. In all other cases,
			 * block_start < block_end should always hold. */
			if (block_start < block_end)
			{
				out_file << prev_ptr->bedoutput(chr_nr, block_start, block_end, graphout) << "\n";
				block_start = block_end + 1;
			}
		}
		
		prev_ptr = ptr;
	}
	
	// report last block
	if (prev_ptr != nullptr)
		out_file << prev_ptr->bedoutput(chr_nr, block_start, prev_ptr->bp_end, graphout) << "\n";
}

/* Converts csv file in_path to a bedgraph file.
 * 
 * If overlap is false, will postprocess overlapping regions
 * to only include the maximum at any given position.
 * */
void csv2bedgraph(const std::string in_path, const std::string out_path, bool overlap)
{
	std::ifstream csv_file;
	csv_file.open(in_path);
	if (!csv_file.is_open())
	{
		std::cerr << "Error opening file " << in_path << std::endl;
		return;
	}
	
	std::string line;
	std::vector<Block> blocks;
	std::getline(csv_file, line); // skip first line
	while (std::getline(csv_file, line))
	{
		Block block = Block(line);
		if (block.sHat <= 0) continue;
		blocks.push_back(block);
	}
	csv_file.close();
	
	std::ofstream out_file;
	out_file.open(out_path, std::ofstream::out | std::ofstream::trunc);
	if (!out_file.is_open())
	{
		std::cerr << "Error opening file " << out_path << std::endl;
		return;
	}
	out_file << "track type=bedGraph name=shat description=selection_coefficient visibility=full graphType=bar windowingFunction=maximum alwaysZero=on\n";
	
	unsigned int chr_nr = guess_chromosome_nr(in_path);
	if (chr_nr == 0)
		std::cerr << "Couldn't get chromosome number from file name " << in_path
			<< ", replace 'chr0' with correct value in output.\n";
	
	if (overlap) // output with possible overlaps
		for (auto block : blocks)
			out_file << "chr" << chr_nr << " " << block << "\n";
	else // remove overlaps
		print_without_overlap(blocks, out_file, chr_nr, true);
	out_file.close();
}


/* Converts csv file in_path to a bed file */
void csv2bed(const std::string in_path, const std::string out_path, bool overlap)
{
	std::ifstream csv_file;
	csv_file.open(in_path);
	if (!csv_file.is_open())
	{
		std::cerr << "Error opening file " << in_path << std::endl;
		return;
	}
	
	std::string line;
	std::vector<Block> blocks;
	std::getline(csv_file, line); // skip first line
	while (std::getline(csv_file, line))
	{
		Block block = Block(line);
		if (block.sHat <= 0) continue;
		blocks.push_back(block);
	}
	csv_file.close();
	
	std::ofstream out_file;
	out_file.open(out_path, std::ofstream::out | std::ofstream::trunc);
	if (!out_file.is_open())
	{
		std::cerr << "Error opening file " << out_path << std::endl;
		return;
	}
	
	unsigned int chr_nr = guess_chromosome_nr(in_path);
	if (chr_nr == 0)
		std::cerr << "Couldn't get chromosome number from file name " << in_path
			<< ", replace 'chr0' with correct value in output.\n";
	
	if (overlap) // output with possible overlaps
		for (auto block : blocks)
			out_file << block.to_bed_line(chr_nr, block.bp_start, block.bp_end) << "\n";
	else // remove overlaps
		print_without_overlap(blocks, out_file, chr_nr);
	out_file.close();
}
