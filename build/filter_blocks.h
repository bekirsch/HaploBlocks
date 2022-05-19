#pragma once
#include <string>

unsigned int filter_blocks(const std::string in_path, const std::string out_path);
void csv2bedgraph(const std::string in_path, const std::string out_path, bool overlap=false);
void csv2bed(const std::string in_path, const std::string out_path, bool overlap=false);
