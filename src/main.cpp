#include <iostream>
#include <string>
//#include "alignment/alignment.hpp"
//#include "minimizers/minimizers.hpp"
//#include "mapper/mapper.hpp"

auto reference = blonde::LoadFasta("reference.fasta");
auto fragments = blonde::LoadFasta("fragments.fasta");
PrintStatistics(reference, fragments);

// Example: alignment
blonde::AlignmentType type = blonde::AlignmentType::Global;
int score = blonde::Align(query_seq, query_len, target_seq, target_len, type, match, mismatch, gap);
std::cout << "Alignment score: " << score << "\n";

// Example: minimizers
auto minimizers = blonde::Minimize(sequence, sequence_len, k, w);
std::cout << "Found " << minimizers.size() << " minimizers\n";

int main(int argc, char** argv) {
    // 1. Parse args (-h, --version, input files, options)
    if (argc < 3) {
        std::cerr << "Usage: blonde_mapper <reference.fasta> <fragments.fasta>\n";
        return 1;
    }

    std::string ref_file = argv[1];
    std::string frag_file = argv[2];

    // 2. Load input files
    auto reference = blonde::LoadFasta(ref_file);
    auto fragments = blonde::LoadFasta(frag_file);

    // 3. Print basic stats
    PrintStatistics(reference, fragments);

    // 4. Alignment test
    auto score = blonde::Align(
        fragments[0].sequence.c_str(), fragments[0].sequence.size(),
        reference[0].sequence.c_str(), reference[0].sequence.size(),
        blonde::AlignmentType::Global,
        1, -1, -1
    );
    std::cout << "Test alignment score: " << score << "\n";

    // 5. Minimizer test
    auto mins = blonde::Minimize(reference[0].sequence.c_str(), reference[0].sequence.size(), 15, 5);
    std::cout << "Reference minimizers: " << mins.size() << "\n";

    return 0;
}