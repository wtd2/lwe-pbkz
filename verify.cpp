#include <NTL/LLL.h>
#include <NTL/ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/matrix.h>
#include <NTL/vec_ZZ.h>

#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

NTL_CLIENT

int main(int argc, char* argv[]) {
    std::string data, sol, fn;
    // std::cin >> data >> sol;
    // std::cin >> fn;
    data = "data/" + std::string(argv[1]) + ".txt";
    sol = "sol/" + std::string(argv[1]) + ".txt";
    ifstream ff, fb;
    ff.open(data.c_str(), ios_base::in);
    mat_ZZ lweA;
    vec_ZZ lweb, lwes;
    int q;
    double sigma2;
    ff >> q;
    ff >> sigma2;
    ff >> lweb;
    ff >> lweA;
    fb.open(sol.c_str(), ios_base::in);
    fb >> lwes;
    double bound = sqrt(lweA.NumRows() * sigma2);
    cout << "lweA: " << lweA.NumRows() << "," << lweA.NumCols() << endl;
    cout << "lweb: " << lweb.length() << endl;
    cout << "lwes: " << lwes.length() << endl;
    mul(lwes, lweA, lwes);
    cout << lwes.length() << endl;

    sub(lweb, lwes, lweb);
    for (int i = 0; i < lweb.length(); i++) {
        lweb[i] = (lweb[i] % q + q) % q;
        if (lweb[i] > q - lweb[i]) lweb[i] -= q;
    }
    double len = 0.0;
    for (int i = 0; i < lweb.length(); i++) len += to_double(lweb[i] * lweb[i]);
    std::cout << lweb << std::endl;
    std::cout << "norm: " << sqrt(len) << std::endl;
    std::cout << "bound: " << bound << std::endl;
    if (sqrt(len) < 3 * bound) {
        return 0;
    }
    return 1;
}
