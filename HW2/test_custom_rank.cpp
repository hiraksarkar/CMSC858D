#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <chrono>
#include "sdsl/bit_vectors.hpp"
#include "custom_bit_op.hpp"

int main(int argc, char** argv){
    if (argc < 3){
        std::cerr << "please provide number of bits and integer \n" ;
        return 1;
    }

    std::cout << "bitarray size " << argv[1] << "\t int value stored " << argv[2] << "\n" ;
    char* end ;
    size_t N = std::strtoul(argv[1], &end, 10) ;
    size_t integer_to_store = std::strtoul(argv[2], &end, 10) ;
    sdsl::bit_vector b(N, 0) ;
    b.set_int(0, integer_to_store, N) ;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    customrank::rank_supp rb(&b) ;
    std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
    std::cout << "Rank Support data structure constructed ... elapsed time " 
               << std::chrono::duration_cast<std::chrono::seconds>(end_t - begin).count() 
               << " sec\n" ;  

    rb.pretty_print() ;
    rb.test_rank() ;
    return 0 ;
}