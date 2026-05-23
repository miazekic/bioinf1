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

long long DiagonalBin(long long diagonal, long long bin_size) {
    if (diagonal >= 0) {
        return diagonal / bin_size;
    }

    return -((-diagonal + bin_size - 1) / bin_size);
}


void test_LIS() {
    std::cout << "\n=== TESTIRANJE LIS ALGORITMA ===" << std::endl;

    
    std::vector<mapper::Hit> hits1 = {
        {10, 500}, {20, 510}, {30, 1200}, {40, 520}, {50, 300}
    };
    auto chain1 = mapper::find_LIS_chain(hits1);
    std::cout << "Test 1 (Buka): ";
    for(auto h : chain1) std::cout << h.target_pos << " ";
    std::cout << (chain1.size() == 3 ? " -> PROŠAO" : " -> PAD") << std::endl;

    
    std::vector<mapper::Hit> hits2 = {
        {1, 100}, {2, 200}, {3, 300}, {4, 400}
    };
    auto chain2 = mapper::find_LIS_chain(hits2);
    std::cout << "Test 2 (Idealno): ";
    for(auto h : chain2) std::cout << h.target_pos << " ";
    std::cout << (chain2.size() == 4 ? " -> PROŠAO" : " -> PAD") << std::endl;

    
    std::vector<mapper::Hit> hits3 = {
        {1, 50}, {2, 40}, {3, 30}, {4, 20}, {5, 10}
    };
    auto chain3 = mapper::find_LIS_chain(hits3);
    std::cout << "Test 3 (Obrnuto): ";
    for(auto h : chain3) std::cout << h.target_pos << " ";
    std::cout << (chain3.size() == 1 ? " -> PROŠAO" : " -> PAD") << std::endl;

    std::cout << "===============================\n" << std::endl;
}


int main(int argc, char** argv) {
    auto t_total_begin = std::chrono::steady_clock::now();

    if (argc < 3) {
    std::cerr << "Usage: ./mapper <reference.fasta> <fragments.fasta>\n";
    return 1;
}

    //auto ref_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>("../../data/referentna.fna");
    //auto sekv = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>("../../data/sekvanca.fasta");
    auto ref_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[1]);
    auto sekv = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[2]);
    
    //std::vector<std::unique_ptr<Sequence>> sequences;
    auto ref = ref_parser->Parse(-1);
    auto sek= sekv->Parse(-1);

        // U main.cpp testnom dijelu:
/*std::string test_query = "GAAAAT";
std::string test_target = "GAT";
std::string cigar;
unsigned int target_begin = 0;

// Obavezno koristi (unsigned int) ispred length() jer length() vraća size_t
int score = mapper::Align(
    test_query.c_str(), (unsigned int)test_query.length(),
    test_target.c_str(), (unsigned int)test_target.length(),
    mapper::AlignmentType::Global,
    2, -1, -1,
    &cigar,
    &target_begin
);

std::cout << "\n--- BRZI TEST ALIGNMENTA ---" << std::endl;
std::cout << "Query:  " << test_query << std::endl;
std::cout << "Target: " << test_target << std::endl;
std::cout << "Score:  " << score << std::endl;
std::cout << "CIGAR:  " << cigar << std::endl;
std::cout << "---------------------------\n" << std::endl;

    std::unordered_map<unsigned int, std::vector<unsigned int>> index2;

    for (const auto& seq : sek) {
        const std::string& s = seq->data;
        auto mins = mapper::Minimize(
            s.c_str(),
            (unsigned int)s.size(),
            5,
            15
        );
        for (const auto& m : mins) {
            unsigned int hash = std::get<0>(m);
            unsigned int pos  = std::get<1>(m);

            index2[hash].push_back(pos);
        }
    }
*/
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
/*
std::cerr << "Reference distinct minimizers: "
          << sorted_frequencies.size() << "\n";

std::cerr << "Ignored frequent minimizers: "
          << remove_count << "\n";

std::cerr << "Indexed distinct minimizers: "
          << index.size() << "\n";
*/

auto t_index_end = std::chrono::steady_clock::now();
/*
std::cerr << "Index build time: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(
                 t_index_end - t_index_begin
             ).count()
          << " ms\n";
*/
for (int i = 0; i < (int)sek.size(); i++) {
    const std::string& frag_data = sek[i]->data;

    // 1. MINIMIZE
    auto frag_mins = mapper::Minimize(
        frag_data.c_str(),
        (unsigned int)frag_data.size(),
        5, 15
    );


    auto t_hits_begin = std::chrono::steady_clock::now();
    // prikupi hitove iz najboljeg dijagonalnog područja
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
    std::cerr << "Fragment " << i
              << " length=" << frag_data.size()
              << " minimizers=" << frag_mins.size()
              << " hits=0\n";
    continue;
}

long long best_bin = 0;
unsigned int best_bin_count = 0;

for (const auto& item : diagonal_counts) {
    if (item.second > best_bin_count) {
        best_bin = item.first;
        best_bin_count = item.second;
    }
}

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

    auto t_hits_end = std::chrono::steady_clock::now();
/*
    std::cerr << "Collect hits time: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(
                 t_hits_end - t_hits_begin
             ).count()
          << " ms\n";

    //debug ispis hits
    std::cerr << "Fragment " << i
          << " length=" << frag_data.size()
          << " minimizers=" << frag_mins.size()
          << " best_diagonal_bin_hits=" << best_bin_count
          << " hits=" << hits.size()
          << "\n";*/


    if (hits.empty()) continue;

    // 2. LIS
    auto t_lis_begin = std::chrono::steady_clock::now();

    auto chain = mapper::find_LIS_chain(hits);

    auto t_lis_end = std::chrono::steady_clock::now();
/*
    std::cerr << "LIS time: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(
                 t_lis_end - t_lis_begin
             ).count()
          << " ms\n";*/

    if (chain.empty()) continue;

    unsigned int q_begin = chain.front().query_pos;
    unsigned int q_end   = chain.back().query_pos + 5;
    unsigned int t_begin = chain.front().target_pos;
    unsigned int t_end   = chain.back().target_pos + 5;
    // osiguraj da regija nije prevelika
    //if (q_end - q_begin > 10000) q_end = q_begin + 10000;
    //if (t_end - t_begin > 10000) t_end = t_begin + 10000;

    // osiguraj granice
    q_end = std::min(q_end, (unsigned int)frag_data.size());
    t_end = std::min(t_end, (unsigned int)ref[0]->data.size());
    // 3. ALIGN
    // TODO obrisi
    //if (q_end <= q_begin || t_end <= t_begin) continue;
    //if (q_end - q_begin > 3000 || t_end - t_begin > 3000) continue;
    std::string cigar;
    bool has_alignment = false;

    if (q_end <= q_begin || t_end <= t_begin) {
        continue;
    }

    if (q_end - q_begin > 12000 || t_end - t_begin > 12000) {
        std::cerr << "Skipping alignment: region too large "
                << "query=" << (q_end - q_begin)
                << " target=" << (t_end - t_begin)
                << "\n";
    } else {
        auto t_align_begin = std::chrono::steady_clock::now();

        mapper::Align(
            frag_data.c_str() + q_begin, q_end - q_begin,
            ref[0]->data.c_str() + t_begin, t_end - t_begin,
            mapper::AlignmentType::Global,
            2, -1, -2,
            &cigar
        );

        auto t_align_end = std::chrono::steady_clock::now();

        std::cerr << "Alignment time: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(
                        t_align_end - t_align_begin
                    ).count()
                << " ms\n";

        has_alignment = true;
    }
    // 4. PAF
    /*{
    std::cout << sek[i]->name      << "\t"
              << frag_data.size()    << "\t"
              << q_begin             << "\t"
              << q_end               << "\t"
              << "+"                 << "\t"
              << ref[0]->name        << "\t"
              << ref[0]->data.size() << "\t"
              << t_begin             << "\t"
              << t_end               << "\t"
              << 255                 << "\n";
    }*/

    {
        unsigned int block_len = std::min(q_end - q_begin, t_end - t_begin);

        std::cout << sek[i]->name        << "\t"
                << frag_data.size()    << "\t"
                << q_begin             << "\t"
                << q_end               << "\t"
                << "+"                 << "\t"
                << ref[0]->name        << "\t"
                << ref[0]->data.size() << "\t"
                << t_begin             << "\t"
                << t_end               << "\t"
                << block_len           << "\t"
                << block_len           << "\t"
                << 255;

        if (has_alignment) {
            std::cout << "\tcg:Z:" << cigar;
        }

        std::cout << "\n";
    }

}
    auto t_total_end = std::chrono::steady_clock::now();
    std::cerr << "Total time: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(
                 t_total_end - t_total_begin
             ).count()
          << " ms\n";
    
    return 0;
}
