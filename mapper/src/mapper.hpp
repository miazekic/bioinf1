#pragma once
#include <vector>
#include <string>

namespace mapper {

struct FastaRecord {
    std::string header;
    std::string sequence;
};

std::vector<FastaRecord> LoadFasta(const std::string& path);

void PrintStatistics(const std::vector<FastaRecord>& reference,
                     const std::vector<FastaRecord>& fragments);

}