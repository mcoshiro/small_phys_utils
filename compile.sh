#!/bin/bash

#TODO move to a makefile like a sane person

#clean
rm lib/distribution_analyzer.o
rm lib/hzg_mem.o
rm lib/spu_utils.o
rm lib/mcmc_sampler.o
rm lib/libSmallPhysUtils.so

#build shared library
g++ `root-config --cflags --glibs` -Iinc -Wall -c src/distribution_analyzer.cpp -o lib/distribution_analyzer.o -fPIC 
g++ `root-config --cflags --glibs` -Iinc -Wall -c src/hzg_mem.cpp -o lib/hzg_mem.o -fPIC 
g++ `root-config --cflags --glibs` -Iinc -Wall -c src/spu_utils.cpp -o lib/spu_utils.o -fPIC 
g++ `root-config --cflags --glibs` -Iinc -Wall -c src/mcmc_sampler.cpp -o lib/mcmc_sampler.o -fPIC 
g++ -shared lib/spu_utils.o lib/distribution_analyzer.o lib/hzg_mem.o lib/mcmc_sampler.o -o lib/libSmallPhysUtils.so

#build executables
#g++ `root-config --cflags --glibs` -Iinc/ -Llib/ -lSmallPhysUtils -Wall src/hzg/hzg_distribution_analyzer.cxx -o bin/hzg_distribution_analyzer 

