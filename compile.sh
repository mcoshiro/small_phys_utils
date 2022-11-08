#!/bin/bash

#TODO move to a makefile like a sane person

#clean
rm lib/distribution_analyzer.o
rm lib/spu_utils.o
rm lib/libSmallPhysUtils.so

#build shared library
g++ `root-config --cflags --glibs` -Iinc -Wall -c src/distribution_analyzer.cpp -o lib/distribution_analyzer.o -fPIC 
g++ `root-config --cflags --glibs` -Iinc -Wall -c src/spu_utils.cpp -o lib/spu_utils.o -fPIC 
g++ -shared lib/spu_utils.o lib/distribution_analyzer.o -o lib/libSmallPhysUtils.so

#build executables
#g++ `root-config --cflags --glibs` -Iinc/ -Llib/ -lSmallPhysUtils -Wall src/hzg/hzg_distribution_analyzer.cxx -o bin/hzg_distribution_analyzer 

