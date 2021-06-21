CC	= g++
SRC	= main.cpp
VRFSRC = verify.cpp
CFLAGS	= -fpermissive -std=c++11 -march=native -msse4 -msse4.1 -mavx --param max-inline-insns-single=5000 -ffast-math -flto -fomit-frame-pointer -fprefetch-loop-arrays -funroll-loops
LDFLAGS	= -lgmp -lmpfr -lgsl -lgslcblas -lntl -fopenmp
INCL	= -I.

all:
	${CC} ${CFLAGS} ${SRC} ${INCL} ${LDFLAGS} -o main

verify:
	${CC} ${CFLAGS} ${VRFSRC} ${INCL} ${LDFLAGS} -o verify


