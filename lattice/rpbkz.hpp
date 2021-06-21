#ifndef _inc_pbkz_hpp
#define _inc_pbkz_hpp

//configures
#define _allow_cachefiles
#define _ibeta_approx_wrapper

//boost library
#include <boost/functional/hash.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/random.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/cpp_int.hpp> 
#include <boost/version.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/functional/hash.hpp>
#include <boost/multiprecision/random.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/math/distributions.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>


using namespace boost::multiprecision;


//system and STL
#include <cstdlib>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>     
#include <algorithm>
#include <functional>        
#include <map>
#include <errno.h>
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include <random>
#include <string>
#include <vector>
#include <map>

using namespace std;

//NTL
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>
#include <NTL/HNF.h>

//gsl
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

NTL_CLIENT

//OpenMP
#include <omp.h>

//To define verbose level (for readability in function call)
#define VL0 0
#define VL1 1
#define VL2 2
#define VL3 3
#define VL4 4
#define VL5 5
#define VL6 6
#define VL7 7

typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<30> > bkzfloat;

typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<10> > float10;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<15> > float15;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<20> > float20;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<25> > float25;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<30> > float30;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<35> > float35;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<40> > float40;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<45> > float45;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50> > float50;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<80> > float80;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100> > float100;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<150> > float150;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<200> > float200;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<300> > float300;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<400> > float400;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<500> > float500;


#define debug_display(X)    {};
        
//Fundemental tools
#include "stringtools.cpp"        
#include "timetools.cpp"        
#include "filetools.cpp"        
#include "memorytools.cpp"
#include "latticetools.cpp"        
#include "misc.cpp"

#include <lattice/gen_uni_mat.cpp>
#include <lattice/pfuncwrapper.cpp>
#include <lattice/enumwrapper.cpp>
#include <lattice/reductionwrapper.cpp>
#include <lattice/pbkzwrapper.cpp>


#endif
