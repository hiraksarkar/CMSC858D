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
#include <fstream>
#include <sstream>
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
                //if((*b)[i] != 0) std::cout << i << "\t" << running_size << "\t" << (*b)[i] << "\n"  ;
                //std::exit(2) ;
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

            //std::cout << "running_size " << running_size << "\n" ;
            maxselect = running_size ;
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

        int select(int pos, bool debug = false){ 

            if(debug) std::cout << "select query " << pos << "\n" ;
            int sup_ind = 0 ; 
            int block_ind = 0;
            int select_answer = 0;
            
            // if pos is more than highest Rank 
            
            // check if this is in the first superblock
            bool inFirstSuperblock = false;
            bool inLastSuperblock = false;

            auto first_superblock_rank = Rs.get_int(1*superblock_word_size , superblock_word_size) ;
            if(pos <= first_superblock_rank or Rs.size() == 1){
                inFirstSuperblock = true ;
            }
            // check if this is in the last superblock
            if(!inFirstSuperblock){
                auto last_superblock_rank = Rs.get_int((Rs.size()-1)*superblock_word_size, superblock_word_size) ;
                if(pos > last_superblock_rank){
                    inLastSuperblock = true ;
                    block_ind = Rs.size() - 1;
                }
            }

            int m_val = 0 ;
            if(!inFirstSuperblock and !inLastSuperblock){
                // Find the first number 
                // greater than or equal to the 
                // search value

                int low = 0, mid = 0 ;
                int high = static_cast<int>(Rs.size()) - 1;
                
                while(low <= high){
                    mid = low + ((high - low) >> 2) ;
                    m_val = Rs.get_int(mid * superblock_word_size, superblock_word_size) ;
                    if(m_val >= pos){
                        high = mid - 1 ;
                    }else{
                        low = mid + 1 ;
                    }
                }
                
                if(low < Rs.size()){
                    sup_ind = low - 1;
                }else{
                    sup_ind = Rs.size() - 1 ;
                }
            }

            if(debug) std::cout << "Sup ind " << sup_ind << "\n" ;
            int offset_pos = pos - Rs.get_int(sup_ind * superblock_word_size, superblock_word_size) ;
            if(debug) std::cout << "offset pos " << offset_pos << "\n" ;
            if(offset_pos == 0){
                //end of prev superblock has the answer
                select_answer = sup_ind * superblock_size ;
                if(debug) std::cout << "answer " << select_answer << "\n" ; 
                //return ;
            }

            // inner block binary search

            int num_block_per_superblock = (superblock_size / block_size) ;
            if(debug) std::cout << "num per superblock " << num_block_per_superblock << "\n" ;
            if(num_block_per_superblock > 1){
                // Find the first number 
                // greater than or equal to the 
                // search value

                int low = 0, mid = 0 ;
                int high = num_block_per_superblock - 1;
                
                while(low <= high){
                    mid = low + ((high - low) >> 2) ;
                    m_val = Rb.get_int((sup_ind * num_block_per_superblock + mid) * block_word_size, block_word_size) ;
                    if(debug) std::cout << "m_val " << m_val << "\n" ;
                    if(m_val >= offset_pos){
                        high = mid - 1 ;
                    }else{
                        low = mid + 1 ;
                    }
                }
                
                if(low < num_block_per_superblock){
                    block_ind = low - 1;
                }else{
                    block_ind = num_block_per_superblock - 1 ;
                }
            }

            if(debug) std::cout << "block ind " << block_ind << "\n" ;
            if(debug) std::cout << "Rb " << Rb << "\n" ;
            int within_block_offset = offset_pos - 
                        Rb.get_int((sup_ind * num_block_per_superblock + block_ind) * block_word_size, block_word_size) ;

            if(debug) std::cout << "within block offset " << within_block_offset << "\n" ;
            select_answer = sup_ind * superblock_size + block_ind * block_size ;
            int step_within_block = 0 ;
            int diff = 0 ;
            while(diff < within_block_offset){
                if((*b)[select_answer + step_within_block]){
                    diff += 1 ;
                }
                step_within_block += 1 ;
            }
            select_answer += step_within_block ;

            sdsl::select_support_mcl<> ssl(b) ;
            
            if(debug) std::cout << "answer " << select_answer << "\n" ; 
            if(debug) std::cout << "sdsl answer " << ssl(pos) << "\n" ;
            return select_answer;
        }
        
        int select_0(int pos, bool debug = false){ 

            if(debug) std::cout << "select query " << pos << "\n" ;
            int sup_ind = 0 ; 
            int block_ind = 0;
            int select_answer = 0;
            
            // if pos is more than highest Rank 
            
            // check if this is in the first superblock
            bool inFirstSuperblock = false;
            bool inLastSuperblock = false;

            auto first_superblock_rank = superblock_size - Rs.get_int(1*superblock_word_size , superblock_word_size) ;
            if(pos <= first_superblock_rank or Rs.size() == 1){
                inFirstSuperblock = true ;
            }
            // check if this is in the last superblock
            if(!inFirstSuperblock){
                auto last_superblock_rank = num_of_superblock * superblock_size - Rs.get_int((Rs.size()-1)*superblock_word_size, superblock_word_size) ;
                if(pos > last_superblock_rank){
                    inLastSuperblock = true ;
                    block_ind = Rs.size() - 1;
                }
            }

            int m_val = 0 ;
            if(!inFirstSuperblock and !inLastSuperblock){
                // Find the first number 
                // greater than or equal to the 
                // search value

                int low = 0, mid = 0 ;
                int high = static_cast<int>(Rs.size()) - 1;
                
                while(low <= high){
                    mid = low + ((high - low) >> 2) ;
                    m_val = mid * superblock_size - Rs.get_int(mid * superblock_word_size, superblock_word_size) ;
                    if(m_val >= pos){
                        high = mid - 1 ;
                    }else{
                        low = mid + 1 ;
                    }
                }
                
                if(low < Rs.size()){
                    sup_ind = low - 1;
                }else{
                    sup_ind = Rs.size() - 1 ;
                }
            }

            if(debug) std::cout << "Sup ind " << sup_ind << "\n" ;
            int offset_pos = pos - (sup_ind * superblock_size - Rs.get_int(sup_ind * superblock_word_size, superblock_word_size)) ;
            if(debug) std::cout << "offset pos " << offset_pos << "\n" ;
            if(offset_pos == 0){
                //end of prev superblock has the answer
                select_answer = sup_ind * superblock_size ;
                if(debug) std::cout << "answer " << select_answer << "\n" ; 
                //return ;
            }

            // inner block binary search

            int num_block_per_superblock = (superblock_size / block_size) ;
            if(debug) std::cout << "num per superblock " << num_block_per_superblock << "\n" ;
            if(num_block_per_superblock > 1){
                // Find the first number 
                // greater than or equal to the 
                // search value

                int low = 0, mid = 0 ;
                int high = num_block_per_superblock - 1;
                
                while(low <= high){
                    mid = low + ((high - low) >> 2) ;
                    m_val = mid * block_size - Rb.get_int((sup_ind * num_block_per_superblock + mid) * block_word_size, block_word_size) ;
                    if(debug) std::cout << "m_val " << m_val << "\n" ;
                    if(m_val >= offset_pos){
                        high = mid - 1 ;
                    }else{
                        low = mid + 1 ;
                    }
                }
                
                if(low < num_block_per_superblock){
                    block_ind = low - 1;
                }else{
                    block_ind = num_block_per_superblock - 1 ;
                }
            }

            if(debug) std::cout << "block ind " << block_ind << "\n" ;
            if(debug) std::cout << "Rb " << Rb << "\n" ;
            int within_block_offset = (offset_pos - ((block_ind * block_size) - Rb.get_int((sup_ind * num_block_per_superblock + block_ind) * block_word_size, block_word_size))) ;

            if(debug) std::cout << "within block offset " << within_block_offset << "\n" ;
            select_answer = sup_ind * superblock_size + block_ind * block_size ;
            int step_within_block = 0 ;
            int diff = 0 ;
            while(diff < within_block_offset){
                if(!(*b)[select_answer + step_within_block]){
                    diff += 1 ;
                }
                step_within_block += 1 ;
            }
            select_answer += step_within_block ;

            sdsl::select_support_mcl<> ssl(b) ;
            
            if(debug) std::cout << "answer " << select_answer << "\n" ; 
            if(debug) std::cout << "sdsl answer " << ssl(pos) << "\n" ;
            return select_answer - 1;
        }

        void test_select(){
            std::cout << "*******************Testing select with sdsl**************\n" ;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            sdsl::select_support_mcl<> ssl(b) ;

            std::cout << "max value for select " << maxselect << "\n" ;
            for(int select_query = 1; select_query <= maxselect ; ++select_query){
                auto computed_index = this->select(select_query) ;
                //sdsl is non inclusive
                //computed_index = computed_index - (*b)[computed_index] ;
                if (computed_index - 1 != ssl(select_query)){
                    //std::cout << select_query << "\t" 
                    //          << computed_index - 1<< "\t"
                    //          << ssl(select_query) << "\n" ;
                    std::cerr << "ranks do not match " 
                              << "computed rank " << computed_index - 1
                              << "\t sdsl rank " << ssl(select_query) 
                              << " BUG !!!! \n" ;
                    std::exit(1) ;
                }
                //std::cout << computed_rank << "\t" << rsb(rank_query) << "\n" ;
                _verbose("\rselect_query passed : %d", select_query);
            }
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << "\nTest completed ... elapsed time " 
                << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() 
                << " sec\n" ;  
                std::cout << "\n*******************Testing Ends***********************\n" ;

        }


        void test_select_0(){
           std::cout << "*******************Testing select with sdsl**************\n" ;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            sdsl::select_support_mcl<> ssl(b) ;

            std::cout << "max value for select " << maxselect << "\n" ;
            for(int select_query = 1; select_query <= maxselect ; ++select_query){
                auto computed_index = this->select(select_query) ;
                //sdsl is non inclusive
                //computed_index = computed_index - (*b)[computed_index] ;
                if (computed_index - 1 != ssl(select_query)){
                    //std::cout << select_query << "\t" 
                    //          << computed_index - 1<< "\t"
                    //          << ssl(select_query) << "\n" ;
                    std::cerr << "ranks do not match " 
                              << "computed rank " << computed_index - 1
                              << "\t sdsl rank " << ssl(select_query) 
                              << " BUG !!!! \n" ;
                    std::exit(1) ;
                }
                //std::cout << computed_rank << "\t" << rsb(rank_query) << "\n" ;
                _verbose("\rselect_query passed : %d", select_query);
            }
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::cout << "\nTest completed ... elapsed time " 
                << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() 
                << " sec\n" ;  
                std::cout << "\n*******************Testing Ends***********************\n" ; 
        }
        
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

        size_t get_rank_0(size_t rank_query){
            if (rank_query >= N){
                std::cerr << "Rank is excedding size of the bit array " << N << "\n" ;
                return 0 ;
            }
            if( rank_query < get_rank(rank_query) ){
                return 0 ;
            }else{
                return( rank_query - get_rank(rank_query) ) + 1 ;
            }

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
            std::cout << "Bit Array size in bits "<< sdsl::size_in_bytes(*b)*8 <<"\n" ;
            std::cout << "Superblock Array size in bits "<< sdsl::size_in_bytes(Rs)*8 <<"\n" ;
            std::cout << "block Array size in bits "<< sdsl::size_in_bytes(Rb)*8 <<"\n" ;
            std::cout << "Inner block array size in bits "<< sdsl::size_in_bytes(Rp)*8 <<"\n" ;
        }

        void overload(){
            std::cout << "Total size in bits " << sdsl::size_in_bytes(Rs)*8 + sdsl::size_in_bytes(Rb)*8 + sdsl::size_in_bytes(Rp)*8 << "\n" ;
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

        size_t maxselect{0} ;
};

class wavelet_tree{
    public:
        wavelet_tree(std::string& outDir_){
            outDir = outDir_ ;
        }

        wavelet_tree(){
            std::string s = "alabar_a_la_alabardadd0000dddfzxsdcxzwesdfggseee_dldls9ldlsdlseddddddddddds" ;
            //std::string s = "0167154263" ;
            std::set<char> alphabetset;
            std::vector<char> alphabets ;

            // int id{0} ;
            for(char &c: s){
                if(alphabetset.find(c) == alphabetset.end()){
                    alphabetset.insert(c) ;
                    alphabets.push_back(c) ;
                }
            }
            std::sort(alphabets.begin(), alphabets.end()) ;
            
            for(size_t i = 0; i < alphabets.size(); ++i){
                alphabetMap[alphabets[i]] = i ;
            }

            txtsize = s.size() ;
            size_t alphabetSize = alphabets.size() ;
            size_t logsize = std::floor(std::log2(alphabetSize-1)) + 1 ;
            std::vector<sdsl::bit_vector> T(s.size()) ;

            for(size_t i = 0 ; i < logsize ; ++i){
                bv_vec.push_back(sdsl::bit_vector(s.size(),0)) ;
            }

            // make a histogram table, can be a vector
            std::vector<int> hist(std::pow(2,logsize),0) ;
            std::vector<size_t> alphabetCoded(s.size(), 0) ;
            for(size_t i = 0; i < s.size(); i++){
                auto c = s[i] ;
                alphabetCoded[i] = alphabetMap[c] ;
                sdsl::bit_vector bt(logsize, 0) ;
                bt.set_int(0, alphabetMap[c]) ;
                T[i] = bt ;
                hist[alphabetMap[c]]++ ;
                bv_vec[0][i] = T[i][logsize - 1];
            }

            //for(size_t i = 0 ; i < hist.size(); ++i){
            //    std::cout << i << "\t" << hist[i] << "\n" ;
            //}
           
            std::cout << "\nConstructed wavelet tree \n" ;

            RangeVec.resize(logsize) ;
            SPosVec.resize(logsize) ;
            SPosVec[0] = {0} ; 
            RangeVec[0] = {0} ;
            for(size_t l = logsize - 1 ; l > 0 ; --l){
                // for this level
                // calculate histogram
                std::vector<size_t> SPos(std::pow(2,l), 0) ;
                std::vector<size_t> range(std::pow(2,l), 0) ;
                std::vector<size_t> reallocatedAlphabet(s.size(),0) ;
                //std::cout << "hist---\n" ;
                for(size_t i = 0 ; i < std::pow(2,l); ++i){
                    hist[i] = hist[2*i] + hist[2*i + 1] ;
                    //std::cout << i << "\t" << hist[i] << "\n" ;
                }
                // use that to calculate SPos
                //std::cout << "\nSPos---\n" ;
                for(size_t i = 1 ; i < std::pow(2,l); ++i){
                    SPos[i] = SPos[i-1] + hist[i - 1] ;
                    //std::cout<< i << "\t" << SPos[i] << "\n" ;
                }
                SPosVec[l] = SPos ;
                size_t bitmask = (std::pow(2,logsize) - 1) - ( std::pow(2, logsize - l) - 1 ) ;
                for(size_t i = 0; i < s.size() ;++i){
                    auto prefix_l = (alphabetCoded[i] >> (logsize - l));
                    auto pos = SPos[prefix_l]++ ;
                    bv_vec[l][pos] = T[i][logsize - l - 1] ;
                    reallocatedAlphabet[pos] = alphabetCoded[i] ;
                }
                // store_
                for(size_t i = 1 ; i < std::pow(2,l) ; ++i){
                    if(i < std::pow(2,l) - 1){
                        range[i] = *std::min_element(reallocatedAlphabet.begin() + SPosVec[l][i], reallocatedAlphabet.begin() + SPosVec[l][i+1]) ;
                    }else{
                        range[i] = *std::min_element(reallocatedAlphabet.begin() + SPosVec[l][i], reallocatedAlphabet.end()) ;
                    }
                }
                RangeVec[l] = range ;

            }

            std::cout << "SPos\n" ;
            for(auto& SPos : SPosVec){
                for(auto& s : SPos){
                    std::cout << s << "\t" ;
                }
                std::cout << "\n" ;
            }

            dump_wavelet_tree() ;
            
            std::vector<rank_supp> bv_vec_r ;
            for(size_t l = 0 ; l < logsize ; ++l){
                bv_vec_r.push_back(rank_supp(&bv_vec[l])) ;
            }
            

            for(size_t test_pos = 0 ; test_pos < s.size() ; ++test_pos){
                for(auto& c : alphabetset){
                    auto fr = get_rank(c, test_pos) ;
                    std::cout << c << "\t" << fr << "\n" ;
                }
            }
            std::cout << s.size() << "\n";
            auto fr = get_rank('x', 0) ;
        }


        void test_load_and_rank(){
            
            std::string s = "alabar_a_la_alabarda" ;
            std::set<char> alphabetset;
            std::vector<char> alphabets ;

            // int id{0} ;
            for(char &c: s){
                if(alphabetset.find(c) == alphabetset.end()){
                    alphabetset.insert(c) ;
                    alphabets.push_back(c) ;
                }
            }
            std::sort(alphabets.begin(), alphabets.end()) ;

            std::cout << "example created \n" ;
            load_wave_tree() ;

            size_t logsize = std::floor(std::log2(alphabetMap.size()-1)) + 1 ;

            std::vector<rank_supp> bv_vec_r ;
            for(size_t l = 0 ; l < logsize ; ++l){
                bv_vec_r.push_back(rank_supp(&bv_vec[l])) ;
            }
            

            for(size_t test_pos = 0 ; test_pos < s.size() ; ++test_pos){
                for(auto& c : alphabetset){
                    auto fr = get_rank(c, test_pos) ;
                }
            }
        }

        void split(const std::string& str, std::vector<std::string>& tokens, const std::string& delim){
            const std::string whiteSpace = " " ;
            size_t prev = 0, pos = 0;
            do
            {
                pos = str.find(delim, prev);
                if (pos == std::string::npos) pos = str.length();
                std::string token = str.substr(prev, pos-prev);
                if (!token.empty() and token != whiteSpace) tokens.push_back(token);
                prev = pos + delim.length();
            }
            while (pos < str.length() && prev < str.length());
        }


        void load_wave_tree(){
            //std::string outDir = "./" ;
            std::cout << "Loading wavelet tree \n" ;

            std::string alphabetMapFile = outDir + "/alphabet_map.txt" ;
            std::string posVecFile = outDir + "/pos_vec.txt" ;

            {
                std::ifstream fileStream(alphabetMapFile.c_str()) ;
                std::string line ;
                // txt size in first line
                std::getline(fileStream, line) ;
                line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                txtsize = std::stoul(line) ;

                while(std::getline(fileStream, line)){
                    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                    std::vector<std::string> tokens ;
                    split(line, tokens, "\t") ;
                        if (tokens.size() > 1){
                            alphabetMap[tokens[0][0]] = std::stoul(tokens[1]) ;
                        }
                }
            }
            std::cout << "Loaded alphabet Map...\n" ;
            std::cout << "Size " << alphabetMap.size() << "\n" ;
            std::cout << "txt size " << txtsize << "\n" ;
            size_t logsize = std::floor(std::log2(alphabetMap.size()-1)) + 1 ;

            {
                std::ifstream fileStream(posVecFile.c_str()) ;
                std::string line ;
                std::vector<std::string> tokens ;

                std::getline(fileStream, line) ;
                split(line, tokens, "\t") ;

                SPosVec.resize(logsize) ;
                RangeVec.resize(logsize) ;

                for(size_t l = 0 ; l < logsize ; ++l){
                    std::vector<size_t> SPos ;
                    SPos.resize(std::pow(2,l)) ;
                    for(size_t i = 0 ; i < std::pow(2,l); ++i){
                        SPos[i] = std::stoul(tokens[(std::pow(2,l)-1) + i]) ;
                    }
                    SPosVec[l] = SPos ;
                }
                tokens.clear() ;
                std::getline(fileStream, line) ;
                split(line, tokens, "\t") ;
                for(size_t l = 0 ; l < logsize ; ++l){
                    std::vector<size_t> range ;
                    range.resize(std::pow(2,l)) ;
                    for(size_t i = 0 ; i < std::pow(2,l); ++i){
                        range[i] = std::stoul(tokens[(std::pow(2,l)-1) + i]) ;
                    }
                    RangeVec[l] = range ;
                }

            }

            std::cout << "Loaded vectors...\n" ;

            {
                bv_vec.resize(logsize) ;
                for(size_t l = 0 ; l < logsize ; ++l){
                    std::string file_l = outDir + "/level_" + std::to_string(l) + ".bin" ;
                    std::ifstream fileStream(file_l.c_str()) ;
                    sdsl::bit_vector b_temp ;
                    sdsl::load(b_temp, fileStream) ;
                    bv_vec[l] = b_temp ;
                }
            }

            std::cout << "loaded \n" ;

        }

        void dump_wavelet_tree(){
            std::string outDir = "./" ;
            std::string alphabetMapFile = outDir + "/alphabet_map.txt" ;
            std::string posVecFile = outDir + "/pos_vec.txt" ;

            {
                std::ofstream fileStream(alphabetMapFile.c_str()) ;
                fileStream << txtsize << "\n" ;
                for(auto& sp : alphabetMap){
                    fileStream << sp.first << "\t" << sp.second << "\n" ;
                }
            }

            {
                std::ofstream fileStream(posVecFile.c_str()) ;
                //fileStream << txtsize << "\n" ;
                //std::cout << "SPos range \n" ;
                for(auto& SPos : SPosVec){
                    for(auto pos : SPos){
                        fileStream << pos << "\t" ;
                    }
                    //std::cout << "" ;
                }
                fileStream << "\n" ;
                
                //std::cout << "Alphabet range \n" ;
                for(auto& range : RangeVec){
                    for(auto r : range){
                        fileStream << r << "\t" ;
                    }
                } 
                fileStream << "\n" ;
            }

            {
                for(size_t l = 0 ; l < bv_vec.size(); ++l){
                    std::string file_l = outDir + "/level_" + std::to_string(l) + ".bin" ;
                    std::ofstream fileStream(file_l.c_str()) ;
                    sdsl::serialize(bv_vec[l] ,fileStream) ;
                }
            }

            std::cout << "dumped\n" ;
        }


        size_t get_rank(char to_search, size_t rank_pos, bool debug = false){
            // build the needed structure
            size_t alphabetSize = alphabetMap.size() ;
            size_t logsize = std::floor(std::log2(alphabetSize-1)) + 1 ;
            std::vector<rank_supp> bv_vec_r ;
            //debug = true ;
            if(debug){
                for(size_t l = 0 ; l < logsize ; ++l){
                    std::cout << bv_vec[l] << "\n" ;
                }
            }


            for(size_t l = 0 ; l < logsize ; ++l){
                bv_vec_r.push_back(rank_supp(&bv_vec[l])) ;
            }
            size_t curr_branch = (alphabetMap[to_search] & static_cast<int>(std::pow(2, logsize - 1)) ) >> (logsize - 1) ;
            size_t curr_rank_query = rank_pos ;
            size_t curr_node_ind = 0 ;
            size_t parent_node_ind = 0 ;

            for(size_t l = 1 ; l < logsize; ++l){
                if(debug) std::cout << "level " << l << "\n" ;

                parent_node_ind = curr_node_ind ;
                if(debug) std::cout << "curr branch " << curr_branch << "\n" ;
                if(debug) std::cout << "curr node ind " << curr_node_ind << "\n" ;
                if(debug) std::cout << "curr_rank_query " << curr_rank_query << "\n" ;
                if(!curr_branch){
                    if(curr_node_ind > 0){
                        curr_node_ind = 2 * curr_node_ind ;
                    }
                    if(curr_node_ind != 0){
                        if(l - 1 > 0){
                            auto start_pos_node = SPosVec[l-1][parent_node_ind] ;
                            if(debug) std::cout << " start_pos_node " << start_pos_node 
                                      << "\t" << " curr_rank_query " << curr_rank_query << "\n" ;
                            curr_rank_query = bv_vec_r[l-1].get_rank_0(start_pos_node + curr_rank_query)  ;
                            if(debug) std::cout << "abs query " << curr_rank_query << "\n" ;
                            if(start_pos_node > 0){
                                if(debug) std::cout << "prev query " << bv_vec_r[l-1].get_rank_0(start_pos_node-1) << "\n" ;
                                curr_rank_query = curr_rank_query - bv_vec_r[l-1].get_rank_0(start_pos_node-1) ;
                                if(debug) std::cout << "curr query " << curr_rank_query << "\n" ;
                            }
                        }else{
                            curr_rank_query = bv_vec_r[l-1].get_rank_0(curr_rank_query) ;
                        }
                    }else{
                        curr_rank_query = bv_vec_r[l-1].get_rank_0(curr_rank_query) ;
                    }
                }else{
                    if(curr_node_ind == 0){
                        curr_node_ind = 1 ;
                    }else{
                        curr_node_ind = 2 * curr_node_ind + 1 ; 
                    }
                    if(curr_node_ind != 0){
                        if(l - 1 > 0){
                            auto start_pos_node = SPosVec[l-1][parent_node_ind] ;
                            if(debug) std::cout << bv_vec[l-1] << "\n" ;
                            if(debug) std::cout << " start_pos_node " << start_pos_node 
                                      << "\t" << " curr_rank_query " << curr_rank_query << "\n" ;
                            curr_rank_query = bv_vec_r[l-1].get_rank(start_pos_node + curr_rank_query) ;
                            if(start_pos_node > 0){
                                if(debug) std::cout << " prev rank " << bv_vec_r[l-1].get_rank(start_pos_node-1) << "\n" ;
                                if(debug) std::cout << " curr rank " << curr_rank_query << "\n" ;
                                curr_rank_query = curr_rank_query - bv_vec_r[l-1].get_rank(start_pos_node-1) ;
                                if(debug) std::cout << "curr_rank_query " << curr_rank_query << "\n" ;
                            }
                        }else{
                            curr_rank_query = bv_vec_r[l-1].get_rank(curr_rank_query) ;
                        }
                    }else{
                        curr_rank_query = bv_vec_r[l-1].get_rank(curr_rank_query) ;
                    }
                }
                if(curr_rank_query > 1)
                    curr_rank_query = curr_rank_query - 1 ;
                curr_branch = (alphabetMap[to_search] &  static_cast<int>(std::pow(2, logsize - l - 1))) >> (logsize - l - 1) ;
            }

            if(debug) std::cout << "curr branch " << curr_branch << "\t" ;
            if(debug) std::cout << "curr_node_ind " << curr_node_ind << "\n" ;
            if(!curr_branch){
                if(curr_node_ind != 0){
                    auto start_pos_node = SPosVec[logsize - 1][curr_node_ind] ;
                    curr_rank_query = bv_vec_r[logsize-1].get_rank_0(start_pos_node + curr_rank_query) - bv_vec_r[logsize-1].get_rank_0(start_pos_node-1);
                }else{
                    curr_rank_query = bv_vec_r[logsize-1].get_rank_0(curr_rank_query) ;
                }
            }else{
                if(curr_node_ind != 0){
                    auto start_pos_node = SPosVec[logsize-1][curr_node_ind] ;
                    curr_rank_query = bv_vec_r[logsize-1].get_rank(start_pos_node + curr_rank_query) - bv_vec_r[logsize-1].get_rank(start_pos_node-1);
                }else{
                    curr_rank_query = bv_vec_r[logsize-1].get_rank(curr_rank_query) ;
                }
            }
            if(debug) std::cout << "rank of " << to_search << "," << rank_pos << "\t" << curr_rank_query << "\n" ;
            return curr_rank_query ;

        }



    private:
        std::map<char, int> alphabetMap ;
        std::vector<sdsl::bit_vector> bv_vec ;
        std::vector<std::vector<size_t>> RangeVec ;
        std::vector<std::vector<size_t>> SPosVec ;

        size_t txtsize ;  
        std::string outDir{"./"} ;


};

} // end namespace custom rank

#endif // CUSTOM_BIT_OP_HPP