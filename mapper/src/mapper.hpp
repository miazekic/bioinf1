#pragma once
#include <vector>
#include <string>

namespace mapper {


    struct Hit {
        unsigned int query_pos;
        unsigned int target_pos;
    };

    
    std::vector<Hit> find_LIS_chain(std::vector<Hit>& hits);

struct FastaRecord {
    std::string header;
    std::string sequence;
};

std::vector<FastaRecord> LoadFasta(const std::string& path);

void PrintStatistics(const std::vector<FastaRecord>& reference,
                     const std::vector<FastaRecord>& fragments);

}