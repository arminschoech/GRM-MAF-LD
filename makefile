
CXX = g++
CXXFLAGS = -O3 -Wall -fopenmp -I/n/app/gsl/2.3/include

grm_maf_ld: main.o geno.o grm.o pheno.o bgen.o indi.o snps.o param.o
	$(CXX) $(CXXFLAGS) -L/n/app/gsl/2.3/lib -L/n/app/intel/2016/mkl/lib/intel64 -o grm_maf_ld main.o geno.o grm.o pheno.o bgen.o indi.o snps.o param.o -lgsl -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -lz -Wl,-rpath,/n/app/intel/2016/mkl/lib/intel64:/n/app/intel/2016/mkl/lib/intel64 -Wl,-rpath,/n/app/gsl/2.3/lib:/n/app/gsl/2.3/lib

main.o: main.cpp geno.h grm.h pheno.h bgen.h indi.h snps.h random.h
	$(CXX) $(CXXFLAGS) -c main.cpp

geno.o: geno.cpp geno.h bgen.h
	$(CXX) $(CXXFLAGS) -c geno.cpp

grm.o: grm.cpp grm.h geno.h indi.h
	$(CXX) $(CXXFLAGS) -c grm.cpp

pheno.o: pheno.cpp pheno.h geno.h indi.h random.h
	$(CXX) $(CXXFLAGS) -c pheno.cpp

bgen.o: bgen.cpp bgen.h indi.h snps.h random.h
	$(CXX) $(CXXFLAGS) -c bgen.cpp

indi.o: indi.cpp indi.h
	$(CXX) $(CXXFLAGS) -c indi.cpp

snps.o: snps.cpp snps.h
	$(CXX) $(CXXFLAGS) -c snps.cpp

param.o: param.cpp param.h
	$(CXX) $(CXXFLAGS) -c param.cpp

