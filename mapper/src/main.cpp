#include <iostream>
#include <string>
#include "../alignment/alignment.hpp"
#include "../minimizers/minimizers.hpp"
#include "mapper.hpp"


int main(int argc, char** argv) {

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





    // 1. Parse args (-h, --version, input files, options)
    if (argc < 3) {
        std::cerr << "Usage: mapper_mapper <reference.fasta> <fragments.fasta>\n";
        return 1;
    }


    
    std::cout << "Reference minimizers: ";

    return 0;
}
