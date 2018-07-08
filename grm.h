#ifndef GRM_H_
#define GRM_H_


#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "geno.h"


class Grm
{
public: 
	Grm(int inN, float inAlpha);
	void addGenotypes(Geno*);
	void divideBySnpNumber();
	
	void writeGrmFiles(std::string fileName, Indi* indiFile);
	void printGrm();
	void printGrm(int maxN);


private:
	int n;
	int m;
	float alpha;
	gsl_matrix_float* pGrm; // pointer to individual x SNP GSL matrix; mean-zero centered minor allele counts

	void writeGrmBinFile(std::string fileName);
	void writeGrmIdFile(std::string fileName, Indi* indiFile);
	void writeGrmNFile(std::string fileName);
};




















#endif
