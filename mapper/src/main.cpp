#include <iostream>
#include <string>
#include "../alignment/alignment.hpp"
#include "../minimizers/minimizers.hpp"
#include "../src/mapper.hpp"
#include "mapper.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include <unordered_map>

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

    auto ref_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>("../../data/referentna.fna");
    auto sekv = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>("../../data/sekvanca.fasta");
    //std::vector<std::unique_ptr<Sequence>> sequences;
    auto ref = ref_parser->Parse(-1);
    auto sek= sekv->Parse(-1);

        // U main.cpp testnom dijelu:
std::string test_query = "GAAAAT";
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

    std::unordered_map<unsigned int, std::vector<unsigned int>> index;
    std::unordered_map<unsigned int, std::vector<unsigned int>> index2;
    for (const auto& seq : ref) {
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

            index[hash].push_back(pos);
        }
    }

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

    test_LIS();

    return 0;
}
