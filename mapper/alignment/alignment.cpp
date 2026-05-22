#include "alignment.hpp"
#include <vector>
#include <string>
#include <algorithm> 

namespace mapper {

typedef struct {
    int cost; /* cost of reaching this cell */
    char parent; /* parent cell */
} cell;

int Align(const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap,
    std::string* cigar,
    unsigned int* target_begin){

    if (cigar == nullptr) {
        std::vector<int> previous(target_len + 1);
        std::vector<int> current(target_len + 1);

        for (unsigned int j = 0; j <= target_len; ++j) {
            previous[j] = static_cast<int>(j) * gap;
        }

        for (unsigned int i = 1; i <= query_len; ++i) {
            current[0] = static_cast<int>(i) * gap;

            for (unsigned int j = 1; j <= target_len; ++j) {
                int diag = previous[j - 1] +
                    (query[i - 1] == target[j - 1] ? match : mismatch);
                int up = previous[j] + gap;
                int left = current[j - 1] + gap;

                current[j] = std::max({diag, up, left});
            }

            std::swap(previous, current);
        }

        if (target_begin != nullptr) {
            *target_begin = 0;
        }

        return previous[target_len];
    }

    std::vector<std::vector<cell>> matrix (query_len + 1, std::vector<cell> (target_len + 1));
    
    matrix[0][0].cost = 0;
    matrix[0][0].parent = 'S';

    for(unsigned int i = 1; i <= query_len; i++){
        matrix[i][0].cost = i*gap;
        matrix[i][0].parent = 'U';
    }

    for(unsigned int j = 1; j <= target_len; j++){
        matrix[0][j].cost = j * gap;
        matrix[0][j].parent = 'L';
    }


    for(unsigned int i = 1; i < query_len + 1; i++){
        for(unsigned int j = 1; j < target_len + 1; j++){
            int diag = matrix[i-1][j-1].cost + (query[i-1] == target[j-1] ? match : mismatch);
            int up = matrix[i-1][j].cost + gap;
            int left = matrix[i][j-1].cost + gap;

            matrix[i][j].cost = std::max({diag, up, left});

            if (matrix[i][j].cost == diag) matrix[i][j].parent = 'D';
            else if (matrix[i][j].cost == up) matrix[i][j].parent = 'U';
            else matrix[i][j].parent = 'L'; 
        }
    }

    if(cigar){
        int i = query_len;
        int j = target_len;
        std::string traceback ="";

        while(i != 0 || j != 0){
            char dir = matrix[i][j].parent;
            if(dir == 'D') {
                traceback += 'M';
                i--;
                j--;
            } else if (dir == 'U'){
                traceback += 'I';
                i--;
            } else if (dir == 'L'){
                traceback += 'D';
                j--;
            }

        }

        std::reverse(traceback.begin(), traceback.end());

        std::string final_cigar = "";
        if(!traceback.empty()){
            int count = 1;
            for (size_t k = 1; k <= traceback.length(); k++){
                if (k<traceback.length() && traceback[k] == traceback[k-1]){
                    count++;
                } else {
                    final_cigar += std::to_string(count) + traceback[k-1];
                    count = 1;
                }

            }
        }
        
        *cigar = final_cigar;  
    }

    if (target_begin != nullptr) {
        *target_begin = 0;
    }

    
    return matrix[query_len][target_len].cost;

    }

}
