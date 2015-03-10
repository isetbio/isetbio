CXX ?= g++
CC ?= gcc
CFLAGS = -Wall -Wconversion -O3 -fPIC
LIBS = blas/blas.a
SHVER = 1
OS = $(shell uname)
#LIBS = -lblas

all: train predict svm-scale

lib: linear.o tron.o blas/blas.a ranksvm.o
	if [ "$(OS)" = "Darwin" ]; then \
		SHARED_LIB_FLAG="-dynamiclib -Wl,-install_name,liblinear.so.$(SHVER)"; \
	else \
		SHARED_LIB_FLAG="-shared -Wl,-soname,liblinear.so.$(SHVER)"; \
	fi; \
	$(CXX) $${SHARED_LIB_FLAG} linear.o tron.o blas/blas.a ranksvm.o -o liblinear.so.$(SHVER)

train: tron.o linear.o train.c blas/blas.a ranksvm.o
	$(CXX) $(CFLAGS) -o train train.c tron.o linear.o ranksvm.o $(LIBS)

predict: tron.o linear.o predict.c blas/blas.a ranksvm.o
	$(CXX) $(CFLAGS) -o predict predict.c tron.o linear.o ranksvm.o $(LIBS)

tron.o: tron.cpp tron.h
	$(CXX) $(CFLAGS) -c -o tron.o tron.cpp

linear.o: linear.cpp linear.h
	$(CXX) $(CFLAGS) -c -o linear.o linear.cpp

blas/blas.a: blas/*.c blas/*.h
	make -C blas OPTFLAGS='$(CFLAGS)' CC='$(CC)';

ranksvm.o: ranksvm.cpp ranksvm.h
	$(CXX) $(CFLAGS) -c -o ranksvm.o ranksvm.cpp

svm-scale:
	$(CXX) $(CFLAGS) -o svm-scale svm-scale.c

clean:
	make -C blas clean
	rm -f *~ tron.o linear.o ranksvm.o train predict liblinear.so.$(SHVER) svm-scale
