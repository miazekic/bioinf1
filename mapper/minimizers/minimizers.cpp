#include "minimizers.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <deque>
namespace mapper{
    static char komplement(char a){
        if(a=='A') return 'T';
        if(a=='T') return 'A';
        if(a=='C') return 'G';
        if(a=='G') return 'C';
        return 'N';
    }

    static char broj(char a){
        if(a=='A') return '0';
        if(a=='T') return '3';
        if(a=='C') return '1';
        if(a=='G') return '2';
        return 0;
    }

static std::string get_complement(std::string sequence, unsigned int len){
    std::string ret(len, ' ');
    for(unsigned int i=0; i<len; i++){
        ret[len-i-1]=komplement(sequence[i]);
    }
    return ret;
}

static std::string get_number(std::string sequence, unsigned int len){
    std::string ret(len, ' ');
    for(unsigned int i=0; i<len; i++){
        ret[i]=broj(sequence[i]);
    }
    return ret;
}

std::vector<std::tuple<unsigned int, unsigned int, bool>> Minimize(
    const char* sequence, unsigned int sequence_len,
    unsigned int kmer_len,
    unsigned int window_len)
{
    
    std::string complement=get_number(get_complement(sequence, sequence_len), sequence_len);
    unsigned n=sequence_len-window_len+1;
    std::string s(sequence);
    std::string sequence2=get_number(s, sequence_len);
    int added=-1;
    std::vector<std::tuple<unsigned int, unsigned int, bool>> result;

    for(unsigned int i=0; i<n; i++){
        
        std::string normal=sequence2.substr(i, window_len);
        //std::string komplementirani(complement.c_str()+sequence_len-window_len-i, window_len);
        bool prva=true;
        std::string min;
        unsigned int j=0;
        unsigned int pos=0;
        for(unsigned int j=0; j<window_len-kmer_len+1; j++){
            if(prva==true){
                min=normal.substr(j, kmer_len);
                prva=false;
                pos=i+j;
            } 
            else{
                if(normal.substr(j, kmer_len)<min){
                    min=normal.substr(j, kmer_len);
                    pos=i+j;
                }
            }
        }
        prva=true;
        if(int(pos)!=added){
            unsigned int h = (unsigned int)std::hash<std::string>{}(min);
            result.emplace_back(h, pos, true);
            added=int(pos);
        }
    }

    return result;
}
}