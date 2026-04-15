#include <iostream>
#include <string>
#include "../alignment/alignment.hpp"
#include "../minimizers/minimizers.hpp"
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

    return 0;
}
