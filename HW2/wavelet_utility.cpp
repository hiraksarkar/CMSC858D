//#include "CLI/CLI.hpp"
#include "clipp.h"

#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <chrono>
#include "sdsl/bit_vectors.hpp"
#include "custom_bit_op.hpp"
#include "ghc/filesystem.hpp"

struct WaveletTreeOpt{
  std::string outDir{""} ;
  std::string referenceFile{""} ;
  std::string rankFile{""} ;
  size_t pos{0} ;
} ;

void buildWaveletTree(WaveletTreeOpt& wtOpt){
    std::string outDir = wtOpt.outDir ;
    if (ghc::filesystem::exists(outDir.c_str())) {
      if (!ghc::filesystem::is_directory(outDir.c_str())) {
        std::cerr << outDir.c_str() << " exists as a file. Cannot create a directory of the same name";
        std::exit(1);
      }
    } else {
      ghc::filesystem::create_directories(outDir.c_str());
    }

    if(!ghc::filesystem::exists(wtOpt.referenceFile)){
        std::cout << "Input file " << wtOpt.referenceFile << " does not exist\n" ; 
    }

    std::ifstream fileStream(wtOpt.referenceFile.c_str()) ;
    std::string line ;
    std::getline(fileStream, line) ;
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());

    customrank::wavelet_tree wt(line) ;
    wt.dump_wavelet_tree(outDir) ;

}

void performRank(WaveletTreeOpt& wtOpt){
    if(!ghc::filesystem::exists(wtOpt.outDir)){
        std::cout << "First call build command to build the index " << wtOpt.outDir << " does not exist\n" ; 
        std::exit(1) ;
    }
    customrank::wavelet_tree wt ;
    wt.load_wave_tree(wtOpt.outDir) ;
    if(!ghc::filesystem::exists(wtOpt.rankFile)){
        std::cout << "Rank file " << wtOpt.rankFile << " does not exist\n" ; 
        std::exit(1) ;
    }
    
    std::ifstream fileStream(wtOpt.rankFile.c_str()) ;
    std::string line ;
    while(std::getline(fileStream , line)){
        std::vector<std::string> tokens ;
        wt.split(line, tokens, "\t") ;
        if(tokens.size() > 1){
            char c = tokens[0][0] ;
            size_t pos = std::stoul(tokens[1]) ;
            if(wt.inAlphabet(c)){
                std::cout<< wt.get_rank(c, pos) << "\n";
            }else{
                std::cout << c << " is not present in the index\n" ;
            }
        }
    }

    //wt.access(21, true) ;

}

void access(WaveletTreeOpt& wtOpt){
    if(!ghc::filesystem::exists(wtOpt.outDir)){
        std::cout << "First call build command to build the index " << wtOpt.outDir << " does not exist\n" ; 
        std::exit(1) ;
    }
    customrank::wavelet_tree wt ;
    wt.load_wave_tree(wtOpt.outDir) ;
    
    if(wtOpt.pos >= wt.size()){
        std::cerr << " accessing outside bound, size() " << wt.size() << "\n" ; 
    }else{
        std::cout << wt.access(wtOpt.pos) << "\n" ;
    }
    //wt.access(21, true) ;

}


int main(int argc, char* argv[]) {
    using namespace clipp ;
    using std::cout ;
    enum class mode {help, build, access, rank, select} ;

    mode selected = mode::help ;
    WaveletTreeOpt wtOpt ;

    auto buildMode = (
     command("build").set(selected, mode::build),

     (option("-r", "--input-file") &
      value("input-file", wtOpt.referenceFile)) %
     "reference file",

    (required("-o", "--out-dir") &
      value("output-dir", wtOpt.outDir)) %
     "output directory"
    ) ;


    auto rankMode = (
     command("rank").set(selected, mode::rank),

     (option("-i", "--input-file") &
      value("input-file", wtOpt.outDir)) %
     "reference file",

    (required("-r", "--rank-file") &
      value("rank-file", wtOpt.rankFile)) %
     "output directory"
    ) ;
    
    auto accessMode = (
     command("access").set(selected, mode::access),

     (option("-i", "--input-file") &
      value("input-file", wtOpt.outDir)) %
     "directory",

    (required("-a", "--access") &
      value("access", wtOpt.pos)) %
     "position"
    ) ;

    auto cli = (
    (buildMode |
     rankMode |
     accessMode | 
        command("--help").set(selected,mode::help) |
        command("-h").set(selected,mode::help) |
        command("help").set(selected,mode::help)
    ),
    
        option("-v", "--version").call([]{std::cout << "version 0.1.0\n\n";}).doc("show version")
    );

    decltype(parse(argc, argv, cli)) res;
    try {
        res = parse(argc, argv, cli);
    } catch (std::exception& e) {
        std::cout << "\n\nParsing command line failed with exception: " << e.what() << "\n";
        std::cout << "\n\n";
        std::cout << make_man_page(cli, "validate") << make_man_page(cli, "extract") ;
        return 1;
    }

    if(res){
    switch(selected){
      case mode::build: buildWaveletTree(wtOpt) ; break ; 
      case mode::rank: performRank(wtOpt) ; break ; 
      case mode::access: access(wtOpt) ; break ; 
    case mode::help: std::cout << make_man_page(cli, "build") << make_man_page(cli, "rank") ; 
    }
    }else{
        cout << usage_lines(cli, "help") << '\n';
    }
    std::exit(0) ;



}
