#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <chrono>
#include "sdsl/bit_vectors.hpp"
#include "custom_bf.hpp"
#include <random>

int main(int argc, char** argv){
    size_t num_of_keys = 1000 ;
    double fpr = 0.01 ;
    std::string s = "something" ;
    std::string ss = "somethink" ;

    bf::bloom_filter cbf(num_of_keys, fpr) ;

    cbf.insert(s) ;
    if(cbf.query(s)){
        std::cout << "SUCC\n" ;
    }else{
        std::cout << "FAIL\n" ;
    }
    if(cbf.query(ss)){
        std::cout << "SUCC\n" ;
    }else{
        std::cout << "FAIL\n" ;
    }

    bf::blocked_bloom_filter bbf(num_of_keys, fpr) ;
    
    bbf.insert(s) ;
    if(bbf.query(s)){
        std::cout << "SUCC\n" ;
    }else{
        std::cout << "FAIL\n" ;
    }
    if(bbf.query(ss)){
        std::cout << "SUCC\n" ;
    }else{
        std::cout << "FAIL\n" ;
    }

    std::cout << "testing storing/loading\n" ;

    cbf.dump_bf("test") ;
    bf::bloom_filter lcbf ;
    std::string folder("test") ;
    lcbf.load_bf(folder) ;
    if(lcbf.query(s)){
        std::cout << "SUCC\n" ;
    }else{
        std::cout << "FAIL\n" ;
    }
    if(lcbf.query(ss)){
        std::cout << "SUCC\n" ;
    }else{
        std::cout << "FAIL\n" ;
    }
    
    bbf.dump_bf("test2") ;
    bf::blocked_bloom_filter lbbf ;
    folder = "test2" ;
    lbbf.load_bf(folder) ;
    if(lcbf.query(s)){
        std::cout << "SUCC\n" ;
    }else{
        std::cout << "FAIL\n" ;
    }
    if(lcbf.query(ss)){
        std::cout << "SUCC\n" ;
    }else{
        std::cout << "FAIL\n" ;
    }
    

    return 0 ;
}