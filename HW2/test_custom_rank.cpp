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
    sdsl::bit_vector b(33, 1) ;
    //b.set_int(0, 140003) ;
    b[10] = 0 ;
    b[11] = 0 ;

    customrank::rank_supp rb(&b) ;
    std::cout << b << "\n" ;
    for(size_t i = 1; i <= 2 ; ++i){
        std::cout << rb.select_0(i)  << "\n";
    }

    //std::cout << "stored integer " << b.get_int(0, N) << "\n";

    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //customrank::rank_supp rb(&b) ;
    //std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
    //std::cout << "Rank Support data structure constructed ... elapsed time " 
    //           << std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() 
    //           << " ms\n" ; 

    std::string outDir = "./" ;
    customrank::wavelet_tree wt(outDir) ; 
    wt.test_load_and_rank() ;
    //customrank::wavelet_tree() ;    
    //rb.overload() ;
    //rb.pretty_print() ;
    //rb.test_rank() ;
    //rb.test_select() ;
    
    
    return 0 ;

}