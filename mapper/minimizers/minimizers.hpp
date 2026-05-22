#pragma once

#include <vector>
#include <tuple>
#include <string>

namespace mapper {

std::vector<std::tuple<unsigned int, unsigned int>> Minimize(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len);

}

