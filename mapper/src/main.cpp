#include <iostream>
#include <string>
#include "../alignment/alignment.hpp"
#include "../minimizers/minimizers.hpp"
#include "../src/mapper.hpp"
#include "mapper.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <unordered_set>
#include <cmath>
#include <vector>
#include <tuple>
#include <cstdlib>

struct Sequence {  
        std::string name;
        std::string data;
    Sequence(const char* name, uint32_t name_len,
                const char* data, uint32_t data_len)
            : name(name, name_len),
            data(data, data_len) {}        
        };

struct Fragment {
    std::string name;
    std::string data;

    
    Fragment(const char* name, uint32_t name_len,
                const char* data, uint32_t data_len)
            : name(name, name_len),
            data(data, data_len) {} 
};

struct MappingResult {
    bool mapped = false;
    char strand = '+';

    unsigned int q_begin = 0;
    unsigned int q_end = 0;

    unsigned int t_begin = 0;
    unsigned int t_end = 0;

    unsigned int chain_length = 0;
    unsigned int best_bin_count = 0;

    long long score = -1;
};

long long DiagonalBin(long long diagonal, long long bin_size) {
    if (diagonal >= 0) {
        return diagonal / bin_size;
    }

    return -((-diagonal + bin_size - 1) / bin_size);
}

long long CandidateScore(
    unsigned int q_begin,
    unsigned int q_end,
    unsigned int chain_length,
    unsigned int best_bin_count)
{
    if (q_end <= q_begin || chain_length == 0) {
        return -1;
    }

    unsigned int query_span = q_end - q_begin;

    return static_cast<long long>(query_span) * 1000
         + static_cast<long long>(chain_length) * 10
         + static_cast<long long>(best_bin_count);
}

static char komplement(char a){
        if(a=='A') return 'T';
        if(a=='T') return 'A';
        if(a=='C') return 'G';
        if(a=='G') return 'C';
        return 'N';
    }

static std::string reverse_complement(const std::string& seq) {
    std::string rc(seq.size(), ' ');
    for (unsigned int i = 0; i < seq.size(); i++) {
        rc[seq.size() - 1 - i] = komplement(seq[i]);
    }
    return rc;
}

MappingResult MapOneStrand(
    const std::string& frag_data,
    char strand,
    const std::unordered_map<unsigned int, std::vector<unsigned int>>& index,
    unsigned int target_len,
    int fragment_index)
{
    MappingResult result;
    result.strand = strand;

    auto frag_mins = mapper::Minimize(
        frag_data.c_str(),
        (unsigned int)frag_data.size(),
        5,
        15
    );

    const long long diagonal_bin_size = 500;
    std::unordered_map<long long, unsigned int> diagonal_counts;

    for (const auto& m : frag_mins) {
        unsigned int hash = std::get<0>(m);
        unsigned int frag_pos = std::get<1>(m);

        auto found = index.find(hash);

        if (found != index.end()) {
            for (unsigned int ref_pos : found->second) {
                long long diagonal =
                    static_cast<long long>(ref_pos) - static_cast<long long>(frag_pos);

                long long bin = DiagonalBin(diagonal, diagonal_bin_size);
                ++diagonal_counts[bin];
            }
        }
    }

    if (diagonal_counts.empty()) {
        std::cerr << "Fragment " << fragment_index
                  << " strand=" << strand
                  << " length=" << frag_data.size()
                  << " minimizers=" << frag_mins.size()
                  << " hits=0\n";
        return result;
    }

    long long best_bin = 0;
    unsigned int best_bin_count = 0;

    for (const auto& item : diagonal_counts) {
        if (item.second > best_bin_count) {
            best_bin = item.first;
            best_bin_count = item.second;
        }
    }

    result.best_bin_count = best_bin_count;

    std::vector<mapper::Hit> hits;
    hits.reserve(best_bin_count * 3);

    for (const auto& m : frag_mins) {
        unsigned int hash = std::get<0>(m);
        unsigned int frag_pos = std::get<1>(m);

        auto found = index.find(hash);

        if (found != index.end()) {
            for (unsigned int ref_pos : found->second) {
                long long diagonal =
                    static_cast<long long>(ref_pos) - static_cast<long long>(frag_pos);

                long long bin = DiagonalBin(diagonal, diagonal_bin_size);

                if (std::llabs(bin - best_bin) <= 1) {
                    hits.push_back({frag_pos, ref_pos});
                }
            }
        }
    }

    std::cerr << "Fragment " << fragment_index
              << " strand=" << strand
              << " length=" << frag_data.size()
              << " minimizers=" << frag_mins.size()
              << " best_diagonal_bin_hits=" << best_bin_count
              << " hits=" << hits.size()
              << "\n";

    if (hits.empty()) {
        return result;
    }

    auto chain = mapper::find_LIS_chain(hits);

    if (chain.empty()) {
        return result;
    }

    result.q_begin = chain.front().query_pos;
    result.q_end = chain.back().query_pos + 5;
    result.t_begin = chain.front().target_pos;
    result.t_end = chain.back().target_pos + 5;

    result.q_end = std::min(result.q_end, (unsigned int)frag_data.size());
    result.t_end = std::min(result.t_end, target_len);

    if (result.q_end <= result.q_begin || result.t_end <= result.t_begin) {
        return result;
    }

    result.chain_length = (unsigned int)chain.size();
    result.score = CandidateScore(
        result.q_begin,
        result.q_end,
        result.chain_length,
        result.best_bin_count
    );

    result.mapped = true;
    return result;
}

int main(int argc, char** argv) {
    auto t_total_begin = std::chrono::steady_clock::now();

    if (argc < 3) {
    std::cerr << "Usage: ./mapper <reference.fasta> <fragments.fasta>\n";
    return 1;
}

    auto ref_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[1]);
    auto sekv = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[2]);
    
    auto ref = ref_parser->Parse(-1);
    auto sek= sekv->Parse(-1);

       
const double f = 0.1;

auto t_index_begin = std::chrono::steady_clock::now();
std::unordered_map<unsigned int, std::vector<unsigned int>> index;

for (const auto& seq : ref) {
    auto ref_mins = mapper::Minimize(
        seq->data.c_str(),
        (unsigned int)seq->data.size(),
        5,
        15
    );

    for (const auto& m : ref_mins) {
        unsigned int hash = std::get<0>(m);
        unsigned int pos = std::get<1>(m);

        index[hash].push_back(pos);
    }
}

std::vector<std::pair<unsigned int, unsigned int>> sorted_frequencies;

for (const auto& item : index) {
    sorted_frequencies.push_back({
        item.first,
        (unsigned int)item.second.size()
    });
}

std::sort(
    sorted_frequencies.begin(),
    sorted_frequencies.end(),
    [](const auto& a, const auto& b) {
        return a.second > b.second;
    }
);

unsigned int remove_count =
    static_cast<unsigned int>(std::ceil(sorted_frequencies.size() * f));

if (f > 0.0 && remove_count == 0 && !sorted_frequencies.empty()) {
    remove_count = 1;
}

for (unsigned int i = 0; i < remove_count && i < sorted_frequencies.size(); ++i) {
    index.erase(sorted_frequencies[i].first);
}

std::cerr << "Reference distinct minimizers: "
          << sorted_frequencies.size() << "\n";

std::cerr << "Ignored frequent minimizers: "
          << remove_count << "\n";

std::cerr << "Indexed distinct minimizers: "
          << index.size() << "\n";


auto t_index_end = std::chrono::steady_clock::now();

std::cerr << "Index build time: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(
                 t_index_end - t_index_begin
             ).count()
          << " ms\n";

for (int i = 0; i < (int)sek.size(); i++) {
    const std::string& original_fragment = sek[i]->data;
    std::string rc_fragment = reverse_complement(original_fragment);

    MappingResult plus = MapOneStrand(
        original_fragment,
        '+',
        index,
        (unsigned int)ref[0]->data.size(),
        i
    );

    MappingResult minus = MapOneStrand(
        rc_fragment,
        '-',
        index,
        (unsigned int)ref[0]->data.size(),
        i
    );

    if (!plus.mapped && !minus.mapped) {
        continue;
    }

    MappingResult best;

    if (plus.mapped && !minus.mapped) {
        best = plus;
    } else if (!plus.mapped && minus.mapped) {
        best = minus;
    } else if (minus.score > plus.score) {
        best = minus;
    } else {
        best = plus;
    }

    std::cerr << "Chosen strand for fragment " << i
              << ": " << best.strand
              << " score=" << best.score
              << " chain_length=" << best.chain_length
              << "\n";

    unsigned int output_q_begin = best.q_begin;
    unsigned int output_q_end = best.q_end;

    if (best.strand == '-') {
        output_q_begin = (unsigned int)original_fragment.size() - best.q_end;
        output_q_end = (unsigned int)original_fragment.size() - best.q_begin;
    }

    unsigned int block_len = std::min(
        output_q_end - output_q_begin,
        best.t_end - best.t_begin
    );

    std::cout << sek[i]->name              << "\t"
              << original_fragment.size()  << "\t"
              << output_q_begin            << "\t"
              << output_q_end              << "\t"
              << best.strand               << "\t"
              << ref[0]->name              << "\t"
              << ref[0]->data.size()       << "\t"
              << best.t_begin              << "\t"
              << best.t_end                << "\t"
              << block_len                 << "\t"
              << block_len                 << "\t"
              << 255                       << "\n";
}

    auto t_total_end = std::chrono::steady_clock::now();

    return 0;
}
