#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <chrono>
#include "sdsl/bit_vectors.hpp"
#include "custom_bf.hpp"
#include <random>

void gen_random_26(std::string& s, const int len) {
    static const char alphanum[] =
        "0123456789"
        "abcdefghijklmnopqrstuvwxyz";

    for (int i = 0; i < len; ++i) {
        s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    //s[len] = 0;
}

std::string gen_random_10(const int len) {
    std::string s(len, '*') ;
    static const char alphanum[] =
        "0123456789";

    for (int i = 0; i < len; ++i) {
        s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    return s ;
    //s[len] = 0;
}

int main(int argc, char** argv){

    // generate 1000 random strings
    std::vector<std::string> examples ;
    examples.reserve(1000) ;
    for(size_t i = 0 ; i < 1000 ; i ++){
        examples.push_back(std::to_string(i)) ;
        //examples.push_back(gen_random_10(10)) ;
    }

    // insert keys from 0 - 500
    // case 1. test with 500 - 1000
    // case 2. test with 250 - 750
    // case 3. test with 0 - 500 

    for(size_t num_of_keys = 1000; num_of_keys < 1000000; num_of_keys += 1000){
        for(double fpr = 0.01 ; fpr < .25 ; fpr = fpr + 0.01){
            double empirical_fpr{0} ;
            bf::bloom_filter bbf(num_of_keys, fpr) ;
            // insert
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            for(size_t i = 0 ; i < 500 ; ++i){
                bbf.insert(examples[i]) ;
            }
            std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
            auto dur_1 = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() ; 
            // query case (1)
            begin = std::chrono::steady_clock::now();
            double fpr_1 = 0 ;
            for(size_t i = 500 ; i < 1000 ; ++i){
                if(bbf.query(examples[i])){
                    fpr_1 += 1 ;
                }
            }
            end_t = std::chrono::steady_clock::now();
            auto dur2 = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() ; 
            
            // query case (2)
            begin = std::chrono::steady_clock::now();
            double fpr_2 = 0 ;
            for(size_t i = 250 ; i < 750 ; ++i){
                if(i > 500 and bbf.query(examples[i])){
                    fpr_2 += 1 ;
                }
            }
            end_t = std::chrono::steady_clock::now();
            auto dur3 = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() ; 
            
            // query case (3)
            begin = std::chrono::steady_clock::now();
            double fpr_3 = 0 ;
            for(size_t i = 0 ; i < 500 ; ++i){
                if(!bbf.query(examples[i])){
                    fpr_3 += 1 ;
                }
            }
            end_t = std::chrono::steady_clock::now();
            auto dur4 = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() ; 

            std::cout << num_of_keys 
                      << "\t" << fpr 
                      << "\t" << fpr_1/500.0 
                      << "\t" << fpr_2/500.0
                      << "\t" << fpr_3/500.0 
                      << "\t" << dur1 
                      << "\t" << dur2
                      << "\t" << dur3
                      << "\t" << dur4
                      << "\n" ;
            
        }
    }

    return 0 ;
}