#include <iostream>
#include <bitset>
#include <string>
#include <climits>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "sdsl/bit_vectors.hpp"
#include "sdsl/int_vector.hpp"

using namespace sdsl ;

template<typename T>
void printVec(std::vector<T>& v){
    for(auto s : v){
        std::cout << s << "\t" ;
    }
    std::cout << "\n" ;
}

void printsdslVec(int_vector<>& v){
    
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
    
    //std::cout << N << "\t" << integer_to_store << "\n" ;

    //
    size_t ws = std::floor(std::log2(N)) ;

    bit_vector* bptr ;
    bptr = new bit_vector(N,0) ;
    std::cout << "bptr size " << bptr->size() << "\n" ;

    bit_vector b(N, 0) ;

    // sdsl rank 
    rank_support_v<> rb(&b) ;

    b.set_int(0, integer_to_store, N) ;
    bptr->set_int(0, integer_to_store, N) ;
    // b[8] = 1;

    // check the size once
    // std::cout << b << "\t" << b.size() << "\n" ;
    
    // naive temporary array
    // try to avoid this

    // create a block size of floor(log(n)/2)
    size_t block_size = static_cast<size_t>(std::floor(std::log2(N) / 2)) ;
    size_t real_block_size = block_size + 1 ;
    
    // create superblock size of b*floor(log(n))
    size_t superblock_size = static_cast<size_t>(block_size * std::floor(std::log2(N))) ;
    size_t num_of_superblock = static_cast<size_t>(std::floor(N / superblock_size)) ;
    // pad by 1 if not multiple
    if( (N % superblock_size) != 0){
        num_of_superblock += 1 ;
    }
    size_t num_of_block = static_cast<size_t>(std::floor(N / block_size)) ;
    // pad by 1 if not multiple
    if(N % block_size != 0){
        num_of_block += 1 ;
    }
    
    // These should ideally be also bit arrays
    int_vector<> Rs(num_of_superblock, 0, superblock_size) ;
    int_vector<> Rb(num_of_block, 0, real_block_size) ;

    //std::cout << Rs.size() << "\t" << Rb.size() << "\n" ;
    std::vector<size_t> Rss(num_of_superblock, 0) ;
    std::vector<size_t> Rbs(num_of_block,0) ;

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
            Rs.set_int(i, running_size, superblock_size);
            std::cout << "i/superblock_size: " << i / superblock_size << "\n" ;
            std::cout << "superblock value " << Rs.get_int(i , superblock_size) << "\n";
            Rss[i / superblock_size] = running_size ;
            inner_sum = 0 ;
        }
        // inside super block
        if(i % block_size == 0){
            auto store_ind = i/block_size ;
            Rb.set_int( store_ind * real_block_size , inner_sum, real_block_size);
            std::cout << Rb.get_int( store_ind * real_block_size , real_block_size ) << "\n"   ;
            Rbs[i / block_size] = inner_sum ;
        }

        inner_sum += b[i] ;
        running_size += b[i] ;
    }

    std::cout << "Rs " << Rs << "\n" ; 
    std::cout << "Rb " << Rb << "\n" ; 

    printVec(Rss) ;
    printVec(Rbs) ;

    // this should be a bit vector
    size_t loglogsize = std::log2(block_size) + 1;
    int_vector<> fixedTypeVec(std::pow(2,block_size) * block_size, 0, loglogsize) ;
    // std::vector<size_t> fixedTypeVec(std::pow(2,block_size) * block_size, 0) ;
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
            //std::cout << b2 << "\t" << i << "\t" << __builtin_popcount(shifted_block) << "\n" ;
            fixedTypeVec.set_int((i * block_size + j)*loglogsize , __builtin_popcount(shifted_block), loglogsize);
        }
    }

    //return 1 ;
    std::cout << "fixedTypeVec " << fixedTypeVec << " loglog " << loglogsize << "\n" ;
    std::cout << b << "\n" ;
    for(size_t query_rank = 0 ; query_rank < N ; ++query_rank){
        auto sup_ind = (query_rank / superblock_size) * superblock_size ;
        auto base_block_ind = (query_rank / block_size) ;
        auto offset = query_rank % block_size ;

        auto content_of_block = b.get_int( base_block_ind * block_size, block_size ) ;

        auto computed_rank = Rs.get_int(sup_ind, superblock_size) + Rb.get_int(base_block_ind * real_block_size, real_block_size) 
                                + fixedTypeVec.get_int( (content_of_block * block_size + offset) * loglogsize , loglogsize) ;

        std::cout << "query rank: " << query_rank
                  << "\t" << "sup: " << sup_ind
                  << "\t" << "sup val: " << Rs.get_int(sup_ind, superblock_size)
                  << "\t" << "blo: " << base_block_ind 
                  << "\t" << "blo val: " << Rb.get_int(base_block_ind * real_block_size, real_block_size)
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

    int_vector<> b3(3, 0, 2) ;
    std::cout << b3.size() << "\t" << b3.width() << "\n" ;

    b3.set_int(0,2) ;
    std::cout << b3 << "\t" << b3.size() << "\t" << b3.width() << "\n" ;
    return 0 ;
    
}