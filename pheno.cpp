#include "pheno.h"


Pheno::Pheno(int inN)
{
	n = inN;
	pPhenotypes = gsl_vector_float_calloc(n);
}


gsl_vector_float* Pheno::simulateEffectVector(int snpN, float causalFrac)
{
	gsl_vector_float* pEffectVector = gsl_vector_float_alloc(snpN);

	for(int i = 0; i < snpN; i++)
	{
		gsl_vector_float_set(pEffectVector, i, (float) gsl_ran_gaussian(ran, 1)); 
		// standard deviation of effect sizes of normalized genotypes is set as 10^-6
		// this does not affect the simulation but taken for numerical conveniences
		if(gsl_rng_uniform(ran) > causalFrac)
			gsl_vector_float_set(pEffectVector, i, 0.0);
	}

	return pEffectVector;
}


void Pheno::addGenEffects(Geno* pGenotypeMatrix, float causalFrac)
{
	gsl_vector_float* pEffectVector = simulateEffectVector(pGenotypeMatrix->getMBuffer(), causalFrac);
	gsl_blas_sgemv(CblasNoTrans, 1, pGenotypeMatrix->getGenMatAddress(), pEffectVector, 1, pPhenotypes);
}






void Pheno::printPhenotypes()
{
	for(int i = 0; i < n; i++)
	{
		std::cout << gsl_vector_float_get(pPhenotypes, i) << std::endl;
	}
	std::cout << std::endl;
}


void Pheno::printPhenotypes(int maxN)
{
	for(int i = 0; i < maxN; i++)
	{
		std::cout << gsl_vector_float_get(pPhenotypes, i) << std::endl;
	}
	std::cout << std::endl;
}

double Pheno::getGeneticVariance() 
{
	double phenoMean = 0;
	double phenoSquaredMean = 0;
	double phenoValue;
	double corrFactor = (double) n / ((double) n + 1.0);

	for(int i = 0; i < n; i++) 
	{
		phenoValue = (double) gsl_vector_float_get(pPhenotypes, i);
		phenoMean += phenoValue / (double) n;
		phenoSquaredMean += pow(phenoValue, 2.0) / (double) n;
	}
	double variance = (phenoSquaredMean - pow(phenoMean, 2.0)) * corrFactor;
	return variance;
}




void Pheno::writePhenoFile(double h2, std::string fileName, Indi* indiFile)
{
	std::ofstream outFile;
	std::string suffix = ".phen";
 	outFile.open(fileName + suffix);

	double environmentalStdDev = sqrt( (1.0 / h2 - 1.0) *  getGeneticVariance() );
	float phenoValue; 

	for(int i = 0; i < n; i++) 
	{
		phenoValue = gsl_vector_float_get(pPhenotypes, i) +  (float) gsl_ran_gaussian(ran, environmentalStdDev);
		outFile << indiFile->getIds()->at(i) << " " << indiFile->getIds()->at(i) << " " << phenoValue << std::endl;
	}

	outFile.close();
}








