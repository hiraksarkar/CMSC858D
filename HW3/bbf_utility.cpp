#include "clipp.h"

#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <chrono>
#include "sdsl/bit_vectors.hpp"
#include "custom_bf.hpp"
#include "ghc/filesystem.hpp"

struct BFOpt{
  std::string outDir{""} ;
  std::string keyFile{""} ;
  std::string inDir{""} ;
  std::string queryFile{""} ;
  size_t numOfKeys{0} ;
  double fpr{0.001} ;
} ;



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

void buildbf(BFOpt& bfOpt){
    bf::bloom_filter cbf(bfOpt.numOfKeys, bfOpt.fpr) ;
    // insert keys
    std::ifstream fileStream(bfOpt.keyFile.c_str()) ;
    std::string line ;
    while(std::getline(fileStream , line)){
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        cbf.insert(line) ;
    }

    cbf.dump_bf(bfOpt.outDir) ;
}

void querybf(BFOpt& bfOpt){
    bf::bloom_filter cbf ;
    cbf.load_bf(bfOpt.inDir) ;
    // insert keys
    std::ifstream fileStream(bfOpt.queryFile.c_str()) ;
    std::string line ;
    while(std::getline(fileStream , line)){
        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        if(cbf.query(line)){
            std::cout << line << ":Y" << "\n" ;
        }else{
            std::cout << line << ":N" << "\n" ;
        }
    }
}


int main(int argc, char* argv[]) {
    using namespace clipp ;
    using std::cout ;
    enum class mode {help, build, query} ;

    mode selected = mode::help ;
    BFOpt bfOpt ;

    auto buildMode = (
     command("build").set(selected, mode::build),

     (required("-k", "--key-file") &
      value("key-file", bfOpt.keyFile)) %
     "file with keys",

    (required("-o", "--out-dir") &
      value("output-dir", bfOpt.outDir)) %
     "output directory",
    
    (required("-n", "--num-keys") &
      value("num-keys", bfOpt.numOfKeys)) %
     "number",

    (required("-f", "--fpr-rate") &
      value("fpr-rate", bfOpt.fpr)) %
     "fpr-rate"

    ) ;


    auto queryMode = (
     command("query").set(selected, mode::query),

     (option("-i", "--input-dir") &
      value("input-file", bfOpt.inDir)) %
     "input directory",

    (required("-q", "--query-file") &
      value("rank-file", bfOpt.queryFile)) %
     "query file"
    ) ;


    auto cli = (
    (buildMode |
     queryMode |
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
      case mode::build: buildbf(bfOpt) ; break ; 
      case mode::query: querybf(bfOpt) ; break ; 
    case mode::help: std::cout << make_man_page(cli, "build") << make_man_page(cli, "query") ; 
    }
    }else{
        cout << usage_lines(cli, "help") << '\n';
    }
    std::exit(0) ;



}