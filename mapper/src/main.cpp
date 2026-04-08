#include <iostream>
#include <string>
#include "../alignment/alignment.hpp"
#include "../minimizers/minimizers.hpp"
#include "mapper.hpp"


int main(int argc, char** argv) {
    // 1. Parse args (-h, --version, input files, options)
    if (argc < 3) {
        std::cerr << "Usage: mapper_mapper <reference.fasta> <fragments.fasta>\n";
        return 1;
    }

    
    std::cout << "Reference minimizers: ";

    return 0;
}
