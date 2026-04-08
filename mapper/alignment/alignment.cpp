#include "alignment.hpp"
#include <vector>
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
    std::string* cigar = nullptr,
    unsigned int* target_begin = nullptr)

    std::vector<std::vector<cell>> matrix(query_len + 1, std::vector<cell>(target_len + 1));
    
    matrix[0][0].score=0;
    matrix[0][0].parent='S';

    for(unsigned int i=0, i<query_len; i++){
        matrix[i][0].cost=i*gap;
        matrix[i][0].parent='U';
    }

    for(unsigned int j=0, j<target_len; i++){
        matrix[0][j].cost=j*gap;
        matrix[0][j].parent='L';
    }

   // int diag=-5;
   // int left=-5;
   // int up=-5;

    for(unsigned int i=1; i<query_len+1; i++){
        for(unsigned int j=1; j<target_len+1; j++){
            int diag=matrix[i-1][j-1].cost+(query[i-1] == target[j-1] ? match : mismatch);
            int up=matrix[i][0].cost+gap;
            int left=matrix[0][j].cost+gap;

            matrix[i][j].cost=std::max(diag, up, left);
            matrix
        }
    }




    return 0;
}

}
