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

    std::string outDir(argv[1]) ;

    std::string rank_benchmark = outDir + "/rank_bench.txt" ;
    std::string select_benchmark = outDir + "/select_bench.txt" ;
    std::string wt_rank_benchmark = outDir + "/wt_rank_bench.txt" ;
    std::string wt_select_benchmark = outDir + "/wt_select_bench.txt" ;


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

    std::ofstream wtSelectStream(wt_select_benchmark.c_str()) ;  
    std::ofstream wtRankStream(wt_rank_benchmark.c_str()) ;  

    for(size_t len = 100000 ; len < 10000000; ++len){

        if(len%100000 == 0){
        
        std::string s(len,'*') ;
        //std::cout << "len " << len << "\n" ;
        gen_random_10(s, len) ;
        customrank::wavelet_tree wt(s) ;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        //customrank::rank_supp rb(&b) ;
        auto f = wt.select(s[len/2],4) ;
        std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
         
        auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() ; 
                
        wtSelectStream << f << "\t" << len << "\t" << wt.getAlphabetSize() << "\t" << dur << "\n" ;
        begin = std::chrono::steady_clock::now();
        f = wt.get_rank(s[len/2],4) ;
        //customrank::rank_supp rb(&b) ;
        end_t = std::chrono::steady_clock::now();
        dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() ; 
        wtRankStream << f << "\t" << len << "\t" << wt.getAlphabetSize() << "\t" << dur << "\n" ;


        }

    }

    std::cout << "10 done \n" ;

    
    for(size_t len = 100000 ; len < 10000000; ++len){

        if(len%100000== 0){
        //std::cout << "len " << len << "\n" ;
        std::string s(len,'*') ;
        gen_random_26(s, len) ;
        customrank::wavelet_tree wt(s) ;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        //customrank::rank_supp rb(&b) ;
        auto f = wt.select(s[len/2],4) ;
        std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
         
        auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() ; 
                
        wtSelectStream << f << "\t" << len << "\t" << wt.getAlphabetSize() << "\t" << dur << "\n" ;
         std::string outDir = "problem" ;
         //wt.dump_wavelet_tree(outDir); 
         //std::exit(1);
        begin = std::chrono::steady_clock::now();
        f = wt.get_rank(s[len/2],4) ;
        //customrank::rank_supp rb(&b) ;
        end_t = std::chrono::steady_clock::now();
        dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() ; 
        wtRankStream << f << "\t" << len << "\t" << wt.getAlphabetSize() << "\t" << dur << "\n" ;
        }

    }
    std::cout << "26 done \n" ;

    for(size_t len = 100000 ; len < 10000000; ++len){

        if(len%100000 == 0){
        std::string s(len,'*') ;
        gen_random_62(s, len) ;
        customrank::wavelet_tree wt(s) ;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        //customrank::rank_supp rb(&b) ;
        auto f = wt.select(s[len/2],4) ;
        std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
         
        auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() ; 
                
        wtSelectStream << f << "\t" << len << "\t" << wt.getAlphabetSize() << "\t" << dur << "\n" ;
        begin = std::chrono::steady_clock::now();
        f = wt.get_rank(s[len/2],4) ;
        //customrank::rank_supp rb(&b) ;
        end_t = std::chrono::steady_clock::now();
        dur = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin).count() ; 
        wtRankStream << f << "\t" << len << "\t" << wt.getAlphabetSize() << "\t" << dur << "\n" ;

        }

    }

    

    std::cout << "62 done \n" ;


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