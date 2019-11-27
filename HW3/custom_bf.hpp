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
#include <cstring>
#include <chrono>
#include <fstream>
#include <sstream>
#include <functional>
#include <array>

#include "xxhash.h"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/int_vector.hpp"
#include "ghc/filesystem.hpp"

#define _verbose(fmt, ...) fprintf(stderr, fmt, __VA_ARGS__)

#define LOG2 std::log(2)
#define DENOM std::pow(std::log(2), 2)
#define CACHESIZE 512

// We use the same principle as
// suggested in double-hashing 
// https://en.wikipedia.org/wiki/Double_hashing
// with some assumptions
// to generate k hashes
size_t stringHasher(
    std::array<uint32_t, 2>& hashPair,
    size_t filterSize, 
    size_t i
){
    return ((hashPair[0] + i * hashPair[1])%filterSize) ;
}

void loadHashPair(std::string& key,
                  std::array<uint32_t, 2>& hashPair
){
    char* pt = const_cast<char*>(key.data());
    size_t val = XXH64(static_cast<void*>(pt), key.size() * sizeof(char), 0);
    std::memcpy(hashPair.data(), &val, sizeof(hashPair));
}

size_t getBlockId(std::string& key, size_t filterSize){
    char* pt = const_cast<char*>(key.data());
    size_t val = XXH64(static_cast<void*>(pt), key.size() * sizeof(char), 0);
    return (val % filterSize) ;
}

namespace bf
{
    class bloom_filter{
        public:
            /*
            N -> number_of_distict keys provided
            by user (we don't have control over this)
            numKeys
            M -> bloom filter size we can control 
            this in order to control a nomilan false 
            positive rate. 
            filterSize
            k -> number of hash function
            numHash
            fpr -> false positive rate

            */
            bloom_filter(size_t& N_, double fpr_){
                fpr = fpr_ ;
                N = N_ ;
                // Use the formula to 
                // to deduce number 
                // of hash functions
                M = (-1) * std::ceil((N * std::log(fpr)) /  DENOM) ;
                k = std::ceil(std::log(M/N) * LOG2) ;

                //std::cout << "filter size: " << M << "\n"
                //          << "num of hashes: " << k << "\n" ;
                b = sdsl::bit_vector(M, 0) ;
            }

            // predetermined filter size
            bloom_filter(size_t M_, size_t N_){
                M = M_ ;
                N = N_ ;
                k = std::ceil(std::log(M/N) * LOG2) ;
                b = sdsl::bit_vector(M, 0) ;
            }

            bloom_filter(sdsl::bit_vector b_ , size_t k_){
                b = b_ ;
                M = b.size() ;
                k = k_ ;
            }
            
            bloom_filter(sdsl::bit_vector b_ , size_t N_, bool load){
                N = N_ ;
                b = b_ ;
                M = b.size() ;
                k = std::ceil(std::log(M/N) * LOG2) ;
            }
            
            // Needed While loading
            bloom_filter(){}

            ~bloom_filter(){}

            size_t access_int(size_t loc, size_t word_size){
                return (b.get_int(loc, word_size)) ;
            }

            void insert(std::string& key){
                std::array<uint32_t, 2> hashPair ;
                loadHashPair(key, hashPair) ;
                for(size_t i = 0; i < k ; ++i){
                    size_t pos = stringHasher(hashPair, M, i) ;
                    b[pos] = 1 ;
                }
            }

            bool query(std::string& key){
                std::array<uint32_t, 2> hashPair ;
                loadHashPair(key, hashPair) ;
                for(size_t i = 0; i < k ; ++i){
                    size_t pos = stringHasher(hashPair, M, i) ;
                    if (!b[pos])
                        return false ;
                }
                return true ;
            }

            void dump_bf(std::string outDir = "./"){

                if (ghc::filesystem::exists(outDir.c_str())) {
                    if (!ghc::filesystem::is_directory(outDir.c_str())) {
                        std::cerr << outDir.c_str() 
                        << " exists as a file. Cannot create a directory of the same name";
                        std::exit(1);
                    }
                    } else {
                    ghc::filesystem::create_directories(outDir.c_str());
                }

                std::string outfile = outDir + "/bf.bin" ;
                std::string metafile = outDir + "/meta.txt" ;

                {
                    std::ofstream fileStream(metafile.c_str()) ;
                    fileStream << k << "\n" ;
                }

                {
                    std::ofstream fileStream(outfile.c_str()) ;
                    sdsl::serialize(b, fileStream) ;
                }
            }

            void load_bf(std::string& inDir){
                std::string outfile = inDir + "/bf.bin" ;
                std::string metafile = inDir + "/meta.txt" ;
                {
                    std::ifstream fileStream(outfile.c_str()) ;
                    sdsl::load(b, fileStream) ;
                }
                {
                    std::ifstream fileStream(metafile.c_str()) ;
                    std::string line ;
                    // txt size in first line
                    std::getline(fileStream, line) ;
                    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                    k = std::stoul(line) ;
                    M = b.size() ;
                }
            }
            
        private:
            // User defined
            size_t N ; // number of keys
            double fpr ; // fpr rate
            // To be calculated
            size_t k ; // number of hash functions
            int M ; // filterSize
            sdsl::bit_vector b ;
    };

    class blocked_bloom_filter{
        public:
            blocked_bloom_filter(size_t N_, double fpr_){
                fpr = fpr_ ;
                N = N_ ;
                // Use the formula to 
                // to deduce number 
                // of hash functions
                M = (-1) * std::ceil((N * std::log(fpr)) /  DENOM) ;

                // Calculate number of blocks
                size_t filterSizePerBlock = CACHESIZE ;
                numBlocks = std::ceil(M/CACHESIZE) ;
                numKeysPerBlock = std::ceil(N / numBlocks) ;
                //std::cout << "number of blocks: " << numBlocks << "\n"
                //          << "block size " << filterSizePerBlock << "\n" 
                //          << "keys per block " << numKeysPerBlock << "\n" ;
                bf_vector.reserve(numBlocks) ;
                for(size_t i = 0 ; i < numBlocks ; ++i){
                    bf_vector.push_back(
                        bloom_filter(filterSizePerBlock, numKeysPerBlock)
                     );
                }
            }

            blocked_bloom_filter(){}

            ~blocked_bloom_filter(){}

            void insert(std::string& key){
                size_t blockId = getBlockId(key, numBlocks) ;
                // insert to that block
                bf_vector[blockId].insert(key) ;
            }

            bool query(std::string& key){
                size_t blockId = getBlockId(key, numBlocks) ;
                // insert to that block
                return bf_vector[blockId].query(key) ;
            }
            
            void dump_bf(std::string outDir = "./"){
                
                if (ghc::filesystem::exists(outDir.c_str())) {
                    if (!ghc::filesystem::is_directory(outDir.c_str())) {
                        std::cerr << outDir.c_str() 
                        << " exists as a file. Cannot create a directory of the same name";
                        std::exit(1);
                    }
                    } else {
                    ghc::filesystem::create_directories(outDir.c_str());
                }
                std::string metafile = outDir + "/meta.txt" ;
                std::string outfile = outDir + "/bf.bin" ;

                {
                    // make large bit vector
                    std::ofstream fileStream(metafile.c_str()) ;
                    fileStream << numKeysPerBlock << "\n" ;
                }

                {
                    // create a bit int array with
                    // word size 
                    size_t word_size = 64 ;
                    size_t num_word_in_cache = std::ceil(CACHESIZE / word_size) ;
                    size_t total_number_of_blocks = numBlocks * num_word_in_cache ;
                    sdsl::int_vector<> out(total_number_of_blocks, 0, word_size) ;
                    // copy the memory block 
                    // ideally this should be done
                    // with memcpy
                    // I am not sure
                    for(size_t i = 0 ; i < numBlocks; ++i){
                        // start index of this block
                        size_t start_index = i * CACHESIZE ;
                        for(size_t j = 0 ; j < num_word_in_cache ; ++j){
                            auto val = bf_vector[i].access_int(j * word_size, word_size) ;
                            out.set_int(start_index + j * word_size, val) ;
                        }
                    }
                    // write down out
                    std::ofstream fileStream(outfile.c_str()) ;
                    sdsl::serialize(out, fileStream) ;

                }
            }

            void load_bf(std::string& inDir){
                
                
                std::string outfile = inDir + "/bf.bin" ;
                std::string metafile = inDir + "/meta.txt" ;

                {
                    std::ifstream fileStream(metafile.c_str()) ;
                    std::string line ;
                    // txt size in first line
                    std::getline(fileStream, line) ;
                    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                    numKeysPerBlock = std::stoul(line) ;
                }
                
                {
                    sdsl::int_vector<> in ;
                    std::ifstream fileStream(outfile.c_str()) ;
                    sdsl::load(in, fileStream) ;

                    // allocate space for each for the 
                    // vector
                    size_t word_size = 64 ;
                    size_t tot_size_in_bits = in.size() ;
                    numBlocks = std::ceil(tot_size_in_bits / CACHESIZE) ;
                    size_t num_word_in_cache = CACHESIZE / word_size ;
                    bf_vector.reserve(numBlocks) ;
                    for(size_t i = 0; i < numBlocks ; ++i){
                        // start index of this block
                        size_t start_index = i * CACHESIZE ;
                        sdsl::bit_vector btemp(CACHESIZE, 0) ;
                        for(size_t j = 0 ; j < num_word_in_cache ; ++j){
                            auto val = in.get_int(start_index + j * word_size, word_size) ;
                            btemp.set_int(j*word_size, val) ;
                        }
                        bf_vector.push_back(
                            bloom_filter(btemp, numKeysPerBlock, true)
                        );
                    }
                }
            }

        private:
            std::vector<bloom_filter> bf_vector ;
            size_t numBlocks ;
            size_t numKeysPerBlock ;
            double fpr ;
            size_t M ;
            size_t N ;
            size_t k ;
    };
}

#endif // CUSTOM_BIT_OP_HPP