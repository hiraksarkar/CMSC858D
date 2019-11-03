#ifndef CUSTOM_BIT_OP_HPP
#define CUSTOM_BIT_OP_HPP

#include <iostream>
#include <bitset>
#include <string>
#include <climits>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include "sdsl/bit_vectors.hpp"
#include "sdsl/int_vector.hpp"

#define _verbose(fmt, ...) fprintf(stderr, fmt, __VA_ARGS__)

namespace customrank
{

class rank_supp{
    public:
        // initialize by an existing bit_vector
        rank_supp(sdsl::bit_vector* b_){
            b = b_ ;
            N = b->size() ;

            // calculate block and super block size
            block_size = static_cast<size_t>(std::floor(std::log2(N) / 2)) ;
            superblock_size = static_cast<size_t>(block_size * std::floor(std::log2(N))) ;

            // calculate the numbers superblocks and blocks
            num_of_block = static_cast<size_t>(std::floor(N / block_size)) ;
            if(N % block_size != 0){
                num_of_block += 1 ;
            }
            num_of_superblock = static_cast<size_t>(std::floor(N / superblock_size)) ;
            if( (N % superblock_size) != 0){
                num_of_superblock += 1 ;
            }

            //calculate the highest size of a superblock_word
            //it might contain a number which is same as 
            //(num_of_superblock - 1)* superblock_size
            superblock_word_size = static_cast<size_t>(std::floor( std::log2(superblock_size * (num_of_superblock)) ))
                                                        + 1 ;

            block_word_size = static_cast<size_t>(std::floor( std::log2(superblock_size - block_size) ))
                                                        + 1 ;

            loglogsize = static_cast<size_t>(std::floor( std::log2(block_size) )) + 1 ;
                                                        + 1 ;
            //now we create support data structure for
            //Rs Rb Rp
            Rs = sdsl::int_vector<>(num_of_superblock, 0, superblock_word_size);
            Rb = sdsl::int_vector<>(num_of_block, 0, block_word_size) ;
            Rp = sdsl::int_vector<>(std::pow(2, block_size) * block_size, 0, loglogsize) ;


            //We would fill up the supporting vectors
            size_t running_size = 0 ;
            size_t prev_superblock_sum = 0 ;
            size_t inner_sum = 0 ;
            for(size_t i = 0; i < N ; ++i){
                if(i % superblock_size == 0){
                    Rs.set_int((i/superblock_size)*superblock_word_size , running_size, superblock_word_size);
                    inner_sum = 0 ;
                }
                // inside super block
                if(i % block_size == 0){
                    auto store_ind = i/block_size ;
                    Rb.set_int( (i/block_size)*block_word_size , inner_sum, block_word_size);
                }

                inner_sum += (*b)[i] ;
                running_size += (*b)[i] ;
            }
            // We fillup the constant table
            for(size_t i = 0 ; i < std::pow(2,block_size) ; i++){
                size_t bit_pattern = 0 ;
                for(size_t j = 0; j < block_size ; j++){
                    bit_pattern = std::pow(2,(j+1)) - 1;
                    auto shifted_block =  i & bit_pattern ;
                    Rp.set_int((i * block_size + j)*loglogsize , __builtin_popcount(shifted_block), loglogsize);
                }
            }
            //This completes the build_phase
        }

        ~rank_supp(){}
        
        size_t get_rank(size_t rank_query){
            auto superblock_ind = (rank_query / superblock_size) ;
            auto block_ind = (rank_query / block_size) ;
            auto offset = rank_query % block_size ;

            auto content_of_superblock = Rs.get_int( superblock_ind * superblock_word_size, superblock_word_size) ;
            auto content_of_block = Rb.get_int( block_ind * block_word_size, block_word_size) ;
            auto content_of_bitarray = b->get_int( block_ind*block_size, block_size) ;
            auto inblock_rank = Rp.get_int( (content_of_bitarray * block_size + offset) * loglogsize, loglogsize) ;

            return (content_of_superblock + content_of_block + inblock_rank) ;

        }

        size_t operator[](size_t rank_query){
            if (rank_query >= N){
                std::cerr << "Rank is excedding size of the bit array " << N << "\n" ;
                return 0 ;
            }
            return get_rank(rank_query) ;
        }

        void test_rank(){

            std::cout << "*******************Testing rank with sdsl**************\n" ;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            sdsl::rank_support_v<> rsb(b) ;
            for(size_t rank_query = 0; rank_query < N ; ++rank_query){
                auto computed_rank = (*this)[rank_query] ;
                //sdsl is non inclusive
                computed_rank = computed_rank - (*b)[rank_query] ;
                if (computed_rank != rsb(rank_query)){
                    std::cerr << "ranks do not match " 
                              << "computed rank " << computed_rank
                              << "\t sdsl rank " << rsb(rank_query) 
                              << " BUG !!!! \n" ;
                    std::exit(1) ;
                }
                //std::cout << computed_rank << "\t" << rsb(rank_query) << "\n" ;
                _verbose("\rrank_query passed : %lu", rank_query);
            }
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << "\nTest completed ... elapsed time " 
                << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() 
                << " sec\n" ;  
                std::cout << "\n*******************Testing Ends***********************\n" ;
        }

        void pretty_print(){
            std::cout << "*****************Rank Datastructure Summary**************\n" ;
            std::cout << "Bit Array size in MB "<< sdsl::size_in_mega_bytes(*b) <<"\n" ;
            std::cout << "Superblock Array size in MB "<< sdsl::size_in_mega_bytes(Rs) <<"\n" ;
            std::cout << "block Array size in MB "<< sdsl::size_in_mega_bytes(Rb) <<"\n" ;
            std::cout << "Inner block array size in MB "<< sdsl::size_in_mega_bytes(Rp) <<"\n" ;
            //std::cout<< "Contents of thr arrays\n" ;
            //std::cout<< "bit array: " << *b << "\n" ;
            //std::cout<< "Rs: " << Rs << "\n" ;
            //std::cout<< "Rb: " << Rb << "\n" ;
            //std::cout<< "Rp: " << Rp << "\n" ; 
            //std::cout << "*************************End****************************\n" ;
        }


    private:
        sdsl::bit_vector* b ;
        sdsl::int_vector<> Rs ;
        sdsl::int_vector<> Rb ;
        sdsl::int_vector<> Rp ;

        size_t N{0} ;
        size_t block_size{0} ;
        size_t superblock_size{0} ;

        size_t block_word_size{0} ;
        size_t superblock_word_size{0} ;
        size_t loglogsize{0} ;

        size_t num_of_block{0} ;
        size_t num_of_superblock{0} ;
};

} // end namespace custom rank

#endif // CUSTOM_BIT_OP_HPP