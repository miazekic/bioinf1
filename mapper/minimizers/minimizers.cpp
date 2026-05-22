#include "minimizers.hpp"

#include <functional>
#include <string>
#include <tuple>
#include <vector>

namespace mapper {

static char komplement(char a) {
    if (a == 'A') return 'T';
    if (a == 'T') return 'A';
    if (a == 'C') return 'G';
    if (a == 'G') return 'C';
    return 'N';
}

static char broj(char a) {
    if (a == 'A') return '0';
    if (a == 'T') return '3';
    if (a == 'C') return '1';
    if (a == 'G') return '2';
    return 0;
}

static std::string get_number(const char* sequence, unsigned int len) {
    std::string ret(len, ' ');

    for (unsigned int i = 0; i < len; ++i) {
        ret[i] = broj(sequence[i]);
    }

    return ret;
}

std::vector<std::tuple<unsigned int, unsigned int>> Minimize(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len)
{
    std::vector<std::tuple<unsigned int, unsigned int>> result;

    if (sequence == nullptr || kmer_len == 0 || window_len == 0 ||
        kmer_len > window_len || sequence_len < window_len) {
        return result;
    }

    const unsigned int n = sequence_len - window_len + 1;
    const unsigned int kmers_in_window = window_len - kmer_len + 1;

    std::string sequence2 = get_number(sequence, sequence_len);

    int added = -1;
    result.reserve(n);

    for (unsigned int i = 0; i < n; ++i) {
        unsigned int pos = i;
        std::string min(sequence2, i, kmer_len);

        for (unsigned int j = 1; j < kmers_in_window; ++j) {
            const unsigned int current_pos = i + j;

            if (sequence2.compare(current_pos, kmer_len, min) < 0) {
                min.assign(sequence2, current_pos, kmer_len);
                pos = current_pos;
            }
        }

        if (static_cast<int>(pos) != added) {
            unsigned int h = static_cast<unsigned int>(
                std::hash<std::string>{}(min)
            );

            result.emplace_back(h, pos);
            added = static_cast<int>(pos);
        }
    }

    return result;
}

}