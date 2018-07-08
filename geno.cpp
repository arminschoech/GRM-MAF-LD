#include "geno.h"


Geno::Geno(int inN, int inM, float inAlpha, float inTau) 
{
	n = inN;
	mBuffer = inM;
	mNextToFill = 0;
	alpha = inAlpha;
	tau = inTau;
	filledWithZeros = 0;

	pGenotypes = gsl_matrix_float_calloc(n, mBuffer);
	pTempVec = gsl_vector_float_alloc(n);

}


void Geno::readInGenVec(Bgen* bgenLink)
{
	if(mNextToFill >= mBuffer) {
      std::cerr << "ERROR: Cannot write to genotype matrix. Genotype matrix is full. " << std::endl;
      exit(1);
    }

	float p = bgenLink->getFrequency();
	if(p>0.5) {
		p = 1 - p;
	}

	// calculate SNP frequency p and normalization constant	
	if(bgenLink->getAlleleCount() <= 0 || bgenLink->getAlleleCount() >= 2*n) 
	{
		gsl_vector_float_set_all(pTempVec, 0.0);
		gsl_matrix_float_set_col(pGenotypes, mNextToFill, pTempVec);
		// if all genotypes are the same (AA or BB), fill in zeros in the 
	} 
	else 
	{
		float normalizationConst = pow(2.0*p*(1-p), alpha/2.0) * pow(calculateLDweight(bgenLink->getLLD()), 0.5);
		if(bgenLink->getFrequency() <= 0.5) {
			// create mean-zero centered and alpha normalized genotype vector at SNP #snpN
			gsl_blas_scopy(bgenLink->genVec, pTempVec);
			gsl_vector_float_add_constant(pTempVec, -2.0*p);
			gsl_blas_sscal(normalizationConst, pTempVec);
		} else {
			gsl_blas_scopy(bgenLink->genVec, pTempVec);
			gsl_blas_sscal(-1.0, pTempVec);
			gsl_vector_float_add_constant(pTempVec, 2.0);
			gsl_vector_float_add_constant(pTempVec, -2.0*p);
			gsl_blas_sscal(normalizationConst, pTempVec);
		}

		// write normalized vector to column #mNextToFill in the genotype matrix
		gsl_matrix_float_set_col(pGenotypes, mNextToFill, pTempVec);
	}
	mNextToFill++;
}


void Geno::fillRest()
{
	gsl_vector_float* zeroVec = gsl_vector_float_calloc(n);
	for( ; mNextToFill < mBuffer; mNextToFill++) {
		gsl_matrix_float_set_col(pGenotypes, mNextToFill, zeroVec);
		filledWithZeros++;
	}
}




void Geno::printGenotypes()
{
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < mBuffer; j++) {
			std::cout << gsl_matrix_float_get(pGenotypes, i, j) << " ";
		}
		std::cout << std::endl;
	}
}


void Geno::printGenotypes(int maxN, int maxM)
{
	for(int i = 0; i < maxN; i++) {
		for(int j = 0; j < maxM; j++) {
			std::cout << gsl_matrix_float_get(pGenotypes, i, j) << " ";
		}
		std::cout << std::endl;
	}
}


float Geno::calculateLDweight(float LLD) 
{
	float value = 1 + LLD * tau;
	if(value >= 0 && value <= 2) {
		return value;
	} else if(value < 0) {
		return 0;
	} else {
		return 2;
	}
}


