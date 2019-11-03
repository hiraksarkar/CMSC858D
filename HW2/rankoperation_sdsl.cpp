#include <iostream>
#include <bitset>
#include <string>
#include <climits>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "sdsl/bit_vectors.hpp"

using namespace sdsl ;

template<typename T>
void printVec(std::vector<T>& v){
    for(auto s : v){
        std::cout << s << "\t" ;
    }
    std::cout << "\n" ;
}

int main(int argc, char** argv)
{

    if (argc < 3){
        std::cerr << "Please provide number of bits and integer \n" ;
        return 1;
    }
    std::cout << argv[0] << "\t" << argv[1] << "\n" ;
    char* end ;
    size_t N = std::strtoul(argv[1], &end, 10) ;
    size_t integer_to_store = std::strtoul(argv[2], &end, 10) ;
    std::cout << N << "\t" << integer_to_store << "\n" ;

    //
    size_t ws = std::floor(std::log2(N)) ;
    bit_vector b(N, 0) ;

    // sdsl rank 
    rank_support_v<> rb(&b) ;

    b.set_int(0, integer_to_store, N) ;
    // b[8] = 1;

    // check the size once
    std::cout << b << "\t" << b.size() << "\n" ;
    
    // naive temporary array
    // try to avoid this

    // create a block size of floor(log(n)/2)
    size_t block_size = static_cast<size_t>(std::floor(std::log2(N) / 2)) ;
    // create superblock size of b*floor(log(n))
    size_t superblock_size = static_cast<size_t>(block_size * std::floor(std::log2(N))) ;

    size_t num_of_superblock = static_cast<size_t>(std::floor(N / superblock_size)) ;
    // pad by 1 if not multiple
    if((N % superblock_size) != 0){
        num_of_superblock += 1 ;
    }
    size_t num_of_block = static_cast<size_t>(std::floor(N / block_size)) ;
    // pad by 1 if not multiple
    if((N % block_size) != 0){
        num_of_block += 1 ;
    }
    
    // These should ideally be also bit arrays
    // int_vector Rs(num_of_superblock, 0, superblock_size)
    // int_vector Rb(num_of_block, 0, block_size)
    std::vector<size_t> Rs(num_of_superblock, 0) ;
    std::vector<size_t> Rb(num_of_block,0) ;


    std::cout << "superblock size: " << superblock_size 
              << "\t" << "number of superblocks: " << num_of_superblock
              << "\t" << "block size: " << block_size 
              << "\t" << "number of blocks: " << num_of_block
              << "\n" ;
    
    // fill the superblock
    // fill the consecutive blocks
    // go over bit array once more
    size_t running_size = 0 ;
    size_t prev_superblock_sum = 0 ;
    //size_t innerblock_index = 0 ;
    size_t inner_sum = 0 ;
    for(size_t i = 0; i < N ; ++i){
        if(i % superblock_size == 0){
            Rs[i / superblock_size] = running_size ;
            inner_sum = 0 ;
        }
        // inside super block
        if(i % block_size == 0){
            Rb[i / block_size] = inner_sum ;
        }
        inner_sum += b[i] ;
        
        running_size += b[i] ;
    }

    printVec(Rs) ;
    printVec(Rb) ;

    // this should be a bit vector
    // int_vector fixed(std::pow(2,block_size) * block_size, 0, std::log2(block_size)) ;
    std::vector<size_t> fixedTypeVec(std::pow(2,block_size) * block_size, 0) ;
    for(size_t i = 0 ; i < std::pow(2,block_size) ; i++){
        size_t bit_pattern = 0 ;
        bit_vector b2(block_size,0) ;
        b2.set_int(0,i, block_size) ;
        for(size_t j = 0; j < block_size ; j++){
            bit_pattern = std::pow(2,(j+1)) - 1;

            // because sdsl is weird !!!!
            auto shifted_block =  i & bit_pattern ;
            //std::cout << "-->" <<  i << "\t" << bit_pattern << "\t" << shifted_block << "\n" ;
            //auto arg = std::bitset<block_size>(i) & std::bitset<block_size>(std::pow(2,j+1)-1) ;
            std::cout << b2 << "\t" << i << "\t" << __builtin_popcount(shifted_block) << "\n" ;
            fixedTypeVec[i * block_size + j] = __builtin_popcount(shifted_block) ;
        }
    }

    for(size_t query_rank = 0 ; query_rank < N ; ++query_rank){
        auto sup_ind = query_rank / superblock_size ;
        auto block_ind = query_rank / block_size ;
        auto offset = query_rank % block_size ;

        auto content_of_block = b.get_int(block_ind*block_size, block_size) ;

        auto computed_rank = Rs[sup_ind] + Rb[block_ind] + fixedTypeVec[content_of_block*block_size + offset] ;

        std::cout << "query rank: " << query_rank
                  << "\t" << "sup: " << sup_ind
                  << "\t" << "blo: " << block_ind 
                  << "\t" << "cont of block: " << content_of_block
                  << "\t" << "computed rank: " << computed_rank 
                  << "\t" << "sdsl rank: " << rb(query_rank) << "\n" ;
    }

    bit_vector b2(6,0) ;
    b2.set_int(0,6,3) ;
    auto val = b2.get_int(0,3) ;
    //bit_vector b3() ;
    std::cout << b2 << "\n" ;
    std::cout << val << "\n" ;
    return 0 ;
    
}