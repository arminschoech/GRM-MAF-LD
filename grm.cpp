#include "grm.h"



Grm::Grm(int inN, float inAlpha)
{
	n = inN;
	m = 0;
	alpha = inAlpha;

	pGrm = gsl_matrix_float_calloc(n, n);
}


void Grm::addGenotypes(Geno* pGenotypeMatrix)
{
	gsl_blas_ssyrk (CblasUpper, CblasNoTrans, 1.0, pGenotypeMatrix->getGenMatAddress(), 1.0, pGrm);
	// adding G * G^T to the GRM, where G is the genotype matrix
	// here I assume that G_ij = (g_ij - 2p_j) / [2p_j(1-p_j)]^alpha/2, 
	// where g_ij is the minor allele count and p_j the frequ.
	m += pGenotypeMatrix->getMBuffer();
	m -= pGenotypeMatrix->filledWithZeros;
}


void Grm::divideBySnpNumber()
{
	gsl_matrix_float_scale(pGrm, 1.0 / (float) m);
	// rescale GRM entries by total number of SNPs included
}



void Grm::writeGrmFiles(std::string fileName, Indi* indiFile)
{
	// create names according GCTA GRM convention
	std::string binSuffix = ".grm.bin";
	std::string NSuffix = ".grm.N.bin";
	std::string idSuffix = ".grm.id";


	std::string binName = fileName + binSuffix;
	std::string NName = fileName + NSuffix;
	std::string idName = fileName + idSuffix;

	// print file names
	// std::cout << binName << std::endl;
	// std::cout << NName << std::endl;
	// std::cout << idName << std::endl;
	
	// create GRMs
	writeGrmBinFile(binName);
	writeGrmNFile(NName);
	writeGrmIdFile(idName, indiFile);
}


void Grm::writeGrmBinFile(std::string fileName)
{
	FILE* outFile;
	outFile = fopen(fileName.c_str(), "wb");
	float buffer;

	// writing out GRM as binary in GCTA order: G_1,1;  G_1,2; G_2,2; G_1.3; G_2.3; G_3,3; ... 
	for(int j = 0; j < n; j++) // loop over the column number of GRM
	{
		for(int i = 0; i <= j; i++) // loop over the row number of GRM
		{	
			buffer = gsl_matrix_float_get(pGrm, i, j);
			fwrite(&buffer, 4, 1, outFile);
		}
	}
	fclose(outFile);
}

void Grm::writeGrmNFile(std::string fileName)
{
	FILE* outFile;
	outFile = fopen(fileName.c_str(), "wb");
	float buffer = m;

	// writing out GRM as binary in GCTA order: G_1,1;  G_1,2; G_2,2; G_1.3; G_2.3; G_3,3; ... 
	for(int j = 0; j < n; j++) // loop over the column number of GRM
	{
		for(int i = 0; i <= j; i++) // loop over the row number of GRM
		{	
			fwrite(&buffer, 4, 1, outFile);
		}
	}
	fclose(outFile);
}

void Grm::writeGrmIdFile(std::string fileName, Indi* indiFile)
{
	std::ofstream outFile;
 	outFile.open(fileName);

	for(int i = 0; i < n; i++) 
	{
		outFile << indiFile->getIds()->at(i) << " " << indiFile->getIds()->at(i) << std::endl;
	}
	outFile.close();
}



void Grm::printGrm()
{
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			std::cout << gsl_matrix_float_get(pGrm, i, j) << " ";
		}
		std::cout << std::endl;
	}
}


void Grm::printGrm(int maxN)
{
	for(int i = 0; i < maxN; i++) {
		for(int j = 0; j < maxN; j++) {
			std::cout << gsl_matrix_float_get(pGrm, i, j) << " ";
		}
		std::cout << std::endl;
	}
}

