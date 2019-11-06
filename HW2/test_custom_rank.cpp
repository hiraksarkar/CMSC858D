#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <chrono>
#include "sdsl/bit_vectors.hpp"
#include "custom_bit_op.hpp"
#include <random>

void gen_random_62(std::string& s, const int len) {
    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";

    for (int i = 0; i < len; ++i) {
        s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    //s[len] = 0;
}

void gen_random_10(std::string& s, const int len) {
    static const char alphanum[] =
        "0123456789";

    for (int i = 0; i < len; ++i) {
        s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    //s[len] = 0;
}

void gen_random_26(std::string& s, const int len) {
    static const char alphanum[] =
        "0123456789"
        "abcdefghijklmnopqrstuvwxyz";

    for (int i = 0; i < len; ++i) {
        s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    //s[len] = 0;
}


int main(int argc, char** argv){

    //char* end ;
    //size_t N = std::strtoul(argv[1], &end, 10) ;

    //// create an extremely conjested array
    //// to see the limit
    //for(size_t bit_size = 10 ; bit_size < 1000; ++bit_size){
    //    sdsl::bit_vector b(bit_size, 1) ;
    //    customrank::rank_supp rb(&b) ;

    //    if(bit_size % 1000 == 0){
    //    size_t bits = rb.overload() ;
    //    auto time_req = rb.test_select() ;
    //     std::cout << bit_size << "\t" << bits << "\t" << time_req << "\n" ;

    //    }

    //}
    for(size_t len = 10000 ; len < 1000000; ++len){
        std::string s(len,'*') ;

        if(len%100000 == 0){
        gen_random_10(s, len) ;
        customrank::wavelet_tree wt(s) ;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        //customrank::rank_supp rb(&b) ;
        auto f = wt.get_rank(s[len/2],len-6) ;
        std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
         
        auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() ; 
                
        std::cout << len << "\t" << wt.getAlphabetSize() << "\t" << dur << "\n" ;

        }

    }
    
    for(size_t len = 10000 ; len < 1000000; ++len){
        std::string s(len,'*') ;

        if(len%100000 == 0){
        gen_random_26(s, len) ;
        customrank::wavelet_tree wt(s) ;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        //customrank::rank_supp rb(&b) ;
        auto f = wt.get_rank(s[len/2],len-6) ;
        std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
         
        auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() ; 
                
        std::cout << len << "\t" << wt.getAlphabetSize() << "\t" << dur << "\n" ;

        }

    }

    for(size_t len = 10000 ; len < 1000000; ++len){
        std::string s(len,'*') ;

        if(len%100000 == 0){
        gen_random_62(s, len) ;
        customrank::wavelet_tree wt(s) ;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        //customrank::rank_supp rb(&b) ;
        auto f = wt.get_rank(s[len/2],len-6) ;
        std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
         
        auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() ; 
                
        std::cout << len << "\t" << wt.getAlphabetSize() << "\t" << dur << "\n" ;

        }

    }

    



    //b.set_int(0, 140003) ;

    //std::cout << b << "\n" ;
    //for(size_t i = 1; i <= 2 ; ++i){
    //    std::cout << rb.select_0(i)  << "\n";
    //}

    //std::cout << "stored integer " << b.get_int(0, N) << "\n";

    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //customrank::rank_supp rb(&b) ;
    //std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
    //std::cout << "Rank Support data structure constructed ... elapsed time " 
    //           << std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() 
    //           << " ms\n" ; 

    //std::string outDir = "./" ;

    //customrank::wavelet_tree wt(outDir) ; 
    //wt.test_load_and_rank() ;
    //customrank::wavelet_tree() ;    
    //rb.pretty_print() ;
    //rb.test_select() ;
    
    
    return 0 ;

}