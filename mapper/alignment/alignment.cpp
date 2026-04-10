#include "alignment.hpp"
#include <vector>
#include <string>
#include <algorithm> 

namespace mapper {

typedef struct {
    int cost; /* cost of reaching this cell */
    char parent; /* parent cell */
} cell;
git
int Align(const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap,
    std::string* cigar = nullptr,
    unsigned int* target_begin = nullptr){

    std::vector<std::vector<cell>> matrix (query_len + 1, std::vector<cell> (target_len + 1));
    
    matrix[0][0].score = 0;
    matrix[0][0].parent = 'S';

    for(unsigned int i = 1, i < query_len; i++){
        matrix[i][0].cost = i*gap;
        matrix[i][0].parent = 'U';
    }

    for(unsigned int j = 1, j < target_len; i++){
        matrix[0][j].cost = j * gap;
        matrix[0][j].parent = 'L';
    }


    for(unsigned int i = 1; i < query_len + 1; i++){
        for(unsigned int j = 1; j < target_len + 1; j++){
            int diag = matrix[i-1][j-1].cost + (query[i-1] == target[j-1] ? match : mismatch);
            int up = matrix[i-1][0].cost + gap;
            int left = matrix[0][j-1].cost + gap;

            matrix[i][j].cost = std::max({diag, up, left});

            if (best == diag) matrix[i][j].parent = 'D';
            else if (best == up) matrix[i][j].parent = 'U';
            else matrix[i][j].parent = 'L'; 
        }
    }

    if(cigar){
        int i = query_len;
        int j = target_len;
        std::string traceback ='';

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

        std::string final_cigar = '';
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
