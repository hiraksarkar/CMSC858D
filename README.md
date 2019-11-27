# CMSC858D
Course Home Work
The code depends on SDSL, once that is installed the utility function for wavelet tree can be compiled with 
```
g++ -g -std=gnu++14 -O0 -g3 -ggdb -pedantic -I ../sdsl-lite/build/include  -L ../sdsl-lite/build/lib -o wt  wavelet_utility.cpp -lsdsl
```

The main code is a single header file `custom_bit_op.hpp`
