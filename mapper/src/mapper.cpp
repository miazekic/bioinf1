#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "../alignment/alignment.hpp"
#include "../minimizers/minimizers.hpp"
#include "mapper.hpp"
#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include <unordered_map>

namespace mapper{
    

std::vector<Hit> find_LIS_chain(std::vector<Hit>& hits) {
    if (hits.empty()) return {};

    // 1. Sortiranje po query poziciji (osigurava linearan poredak u fragmentu)
    std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
    if (a.query_pos != b.query_pos) {
        return a.query_pos < b.query_pos;
    }
    return a.target_pos > b.target_pos;
    });

    int n = hits.size();
    std::vector<int> parent(n, -1);
    std::vector<int> tails_indices; // Čuva indekse iz 'hits' vektora
    tails_indices.reserve(n);


    for (int i = 0; i < n; i++) {
        // Binarna pretraga nad target_pos vrijednostima
        auto it = std::lower_bound(tails_indices.begin(), tails_indices.end(), i,
            [&hits](int a_idx, int i_idx) {
                return hits[a_idx].target_pos < hits[i_idx].target_pos;
            });

        int dist = std::distance(tails_indices.begin(), it);

        if (it == tails_indices.end()) {
            if (!tails_indices.empty()) {
                parent[i] = tails_indices.back();
            }
            tails_indices.push_back(i);
        } else {
            if (dist > 0) {
                parent[i] = tails_indices[dist - 1];
            }
            *it = i;
        }
    }

    // 4. Rekonstrukcija (Traceback kroz roditelje)
    std::vector<Hit> chain;
    int curr = tails_indices.back();
    while (curr != -1) {
        chain.push_back(hits[curr]);
        curr = parent[curr];
    }
    std::reverse(chain.begin(), chain.end());

    return chain;
}
}