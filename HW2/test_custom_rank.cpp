#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <chrono>
#include "sdsl/bit_vectors.hpp"
#include "custom_bit_op.hpp"

int main(int argc, char** argv){

    //char* end ;
    //size_t N = std::strtoul(argv[1], &end, 10) ;

    //// create an extremely conjested array
    //// to see the limit
    for(size_t bit_size = 1000 ; bit_size < 10000000 ; ++bit_size){
        sdsl::bit_vector b(bit_size, 1) ;
        customrank::rank_supp rb(&b) ;

        size_t bits = rb.overload() ;
        auto time_req = rb.test_rank() ;
        if(bit_size % 100000 == 0)
            std::cout << bit_size << "\t" << bits << "\t" << time_req << "\n" ;

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