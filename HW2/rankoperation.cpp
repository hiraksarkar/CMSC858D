#include <iostream>
#include <bitset>
#include <string>
#include <climits>
#include <cmath>
#include <vector>
#include <cstdio>

//using namespace std;
template<typename T>
void printVec(std::vector<T>& v){
    for(auto s : v){
        std::cout << s << "\t" ;
    }
    std::cout << "\n" ;
}


int main()
{
    const size_t N = 16 ;
    std::bitset<N> bs(140003) ;
    std::vector<size_t> naiveArray(bs.size(),0) ;
    
    size_t running_sum{0} ;
    for(size_t i = 0 ; i < bs.size() ; ++i){
        
        running_sum += bs[i] ;
        naiveArray[i] = running_sum ;
    }
    
    printVec<size_t>(naiveArray) ;
    
    std::cout << bs << "\n" ;
    
    // block size b = std::floor(log(n)/2) 
    const size_t block_size = static_cast<size_t>(std::floor(std::log2(N) / 2)) ;
    const size_t superblock_size = static_cast<size_t>(block_size * std::floor(std::log2(N))) ;
    
    std::cout << block_size << "\t" << superblock_size << "\n" ;
    
    size_t num_of_superblock = static_cast<size_t>(std::floor(N / superblock_size)) ;
    
    std::vector<size_t> superblockArray(num_of_superblock, 0) ;
    for(size_t i = 0 ; i < num_of_superblock; ++i){
        superblockArray[i] = naiveArray[i * superblock_size] ; 
    }
    
    std::cout << "Super Block \n" ; 
    printVec<size_t>(superblockArray) ;
    
    // crate blockarray
    size_t num_of_block = static_cast<size_t>(std::floor(N / block_size)) ;
    std::vector<size_t> blockArray(num_of_block, 0) ;
    for(size_t k = 0 ; k < num_of_block; k++){
        // which superblock does it belong to
        size_t j = static_cast<size_t>(k / std::floor(std::log2(N))) ;
        blockArray[k] = naiveArray[k * block_size] - naiveArray[j * superblock_size] ;
        // std::cout << naiveArray[k * block_size] << "\t" << blockArray[k] << "\t" << j << "\n" ;
    }
    
    // Create fixed type of size b
    std::cout << "Block \n" ;
    printVec<size_t>(blockArray) ;
    
    std::vector<size_t> fixedTypeVec(std::pow(2,block_size) * block_size, 0) ;
    
    std::cout << "--\n" ; 
    for(size_t i = 0 ; i < std::pow(2,block_size) ; i++){
        for(size_t j = 0; j < block_size ; j++){
            auto arg = std::bitset<block_size>(i) & std::bitset<block_size>(std::pow(2,j+1)-1) ;
            std::cout << i << "\t" << __builtin_popcount(arg.to_ulong()) << "\n" ;
            fixedTypeVec[i * block_size + j] = __builtin_popcount(arg.to_ulong()) ;
        }
    }
    
    printVec<size_t>(fixedTypeVec) ;

    size_t query_rank = 10 ;
    auto sup_r = query_rank / superblock_size ;
    auto block_r = query_rank / block_size ;
    auto offset = query_rank % block_size ;
    
    std::bitset<N> blockmask(std::pow(2,block_size) - 1) ;
    auto shifted = (blockmask << (block_r * block_size)) ;
    auto val = (bs & shifted) ;
    
    std::cout << superblockArray[sup_r] + blockArray[block_r] +  fixedTypeVec[val.to_ulong() * block_size + offset] << "\n" ;
    std::cout << naiveArray[query_rank] << "\n" ;
    // std::bitset<N> blockmask(std::pow(2,block_r+1) - std::pow(2,block_r))
    // std::cout << std::bitset<N>(bs.begin() + block_r, bs.begin() + block_size) << "\n";

    return 0;
}