#pragma once
#include <vector>
#include <string>

namespace mapper {

struct Minimizer {
    size_t position;
    std::string kmer;
};

std::vector<Minimizer> Minimize(const char* sequence, size_t len, int k, int w);

}