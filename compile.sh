#!/bin/bash

#TODO move to a makefile like a sane person

#clean
rm lib/distribution_analyzer.o
rm lib/hzg_mem.o
rm lib/spu_utils.o
rm lib/mcmc_sampler.o
rm lib/correction.o
rm lib/correction_wrapper.o
rm lib/libSmallPhysUtils.so

#build shared library
g++ -o lib/distribution_analyzer.o -c -O2 -isystem `root-config --incdir` `root-config --cflags` -pedantic -Wall -fPIC -Iinc src/distribution_analyzer.cpp
g++ -o lib/hzg_mem.o -c -O2 -isystem `root-config --incdir` `root-config --cflags` -pedantic -Wall -fPIC -Iinc src/hzg_mem.cpp
g++ -o lib/spu_utils.o -c -O2 -isystem `root-config --incdir` `root-config --cflags` -pedantic -Wall -fPIC -Iinc src/spu_utils.cpp
g++ -o lib/mcmc_sampler.o -c -O2 -isystem `root-config --incdir` `root-config --cflags` -pedantic -Wall -fPIC -Iinc src/mcmc_sampler.cpp
g++ -o lib/correction.o -c -O2 -isystem `root-config --incdir` `root-config --cflags` -std=c++17 -pedantic -Wall -fPIC -Iinc src/correction.cpp
g++ -o lib/correction_wrapper.o -c -O2 -isystem `root-config --incdir` `root-config --cflags` -std=c++17 -pedantic -Wall -fPIC -Iinc src/correction_wrapper.cpp
g++ -o lib/libSmallPhysUtils.so `root-config --glibs` `root-config --ldflags` -shared lib/spu_utils.o lib/distribution_analyzer.o lib/hzg_mem.o lib/mcmc_sampler.o lib/correction.o lib/correction_wrapper.o

#build executables
#g++ `root-config --cflags --glibs` -Iinc/ -Llib/ -lSmallPhysUtils -Wall src/hzg/hzg_distribution_analyzer.cxx -o bin/hzg_distribution_analyzer 
#g++ -o bin/hzg_distribution_analyzer -c -O2 -isystem `root-config --incdir` `root-config --cflags` -pedantic -Wall -fPIC -Iinc src/hzg/hzg_distribution_analyzer.cxx

