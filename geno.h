#ifndef GENO_H_
#define GENO_H_


#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "bgen.h"


class Geno
{

public: 
	Geno(int inN, int inM, float inAlpha, float inTau);
	int getMBuffer();

	void resetGenotypeMatrix(); // genotype matrix resets to start filling it from first column onwards
	void readInGenVec(Bgen*); // read in the genotype vector from Bgen object into the next column of the genotype matrix
	void fillRest();
	gsl_matrix_float* getGenMatAddress();
	
	void printGenotypes();
	void printGenotypes(int maxN, int maxM); 
	int filledWithZeros; // number of SNPs filled with zero genotypes due to lack of further SNPs to read in

private:
	int n;
	int mBuffer; // maximal SNP number, i.e. size of the genotype matrix
	int mNextToFill; // column in genotype matrix that is next to fill
	float alpha;
	float tau;
	gsl_matrix_float* pGenotypes; // pointer to individual x SNP GSL matrix; mean-zero centered and alpha normalized
	gsl_vector_float* pTempVec; // temporary vector to store genotypes for current SNP

	float calculateLDweight(float LLD);
};

inline gsl_matrix_float* Geno::getGenMatAddress()
{
	return pGenotypes;
}

inline int Geno::getMBuffer() 
{
	return mBuffer;
}

inline void Geno::resetGenotypeMatrix()
{
	mNextToFill = 0;
}












#endif