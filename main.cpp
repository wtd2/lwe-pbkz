#include "lattice/rpbkz.hpp"
#include "lattice/bkztest.cpp"


int main(int argc, char** argv) {
    if (strcmp(argv[1], "pbkz") == 0) {
        pbkz_lwe(argv[2]);
    } else if (strcmp(argv[1], "kannan") == 0) {
        kannan_lwe(argv[2]);
    } else if (strcmp(argv[1], "deep") == 0) {
        // deep_lwe(argv[2]);
    }
    
    exit(0);
}
