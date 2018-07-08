#ifndef BGEN_H_
#define BGEN_H_


#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <string>
#include <gsl/gsl_vector.h>
#include <zlib.h>

#include "random.h"
#include "indi.h"
#include "snps.h"



class Bgen 
{
public:
	
	Bgen(const char* fileName, int partNum, int totalNum, int inMinimalAlleleCount, Indi* inIndividuals, Snps* inSnpFile);
	// reads in the nth part when splitting all SNPs in the BGEN file into a total of ofTotal parts
	// takes individuals from individuals file

	int getN();
	int getAlleleCount();
	float getFrequency();
	float getLLD();	

	void sampleNextNonZeroSnp();
	void mapNextNonZeroSnp();
	void dosageNextNonZeroSnp();
	bool loadSnpSuccess; // indicates if loading next SNP was successful
	
	gsl_vector_float* genVec; // sampled or MAP genotype vector
	bool reachedEnd; // bool that indicates TRUE if there are no more SNPs to read in
	
	void printGenProb();
	void printGenProb(int maxN);
	void printGenVec(int maxN);


private:
	FILE* bgenFile;
	int nBgen; // number of individuals in the bgen file
	int mBgen; // the total number of SNPs in the BGEN file

	Indi* individuals;
	int n; // number of individuals used (as set by indi object)

	Snps* snpFile;

	int startM; // first SNP to be read in in this simulation
	int endM; // last SNP to be read in in this simulation
	int currM; // current SNP


	unsigned char* zBuf; // buffer for compressed data
  	unsigned short* shortBuf; // buffer for unpacked data

  	int alleleCount; // number of B alleles in genVec
  	int minAlleleCount; // SNPs with lower allele counts are not used
  	float frequency; // frequency of B alleles in genVec
  	float LLD;
  	std::string currRsId;
  	std::string currSnpId;

  	int readSnpHeader(); // reads in SNP header and returns size of compressed SNP data; used by loadNextSnp()
	void loadNextSnp();
	bool useSnp();
	void sampleSnp();
	void mapSnp();
	void dosageSnp();
	void skipSnps(int skipM);

};


inline int Bgen::getN()
{
	return n;
}

inline int Bgen::getAlleleCount()
{
	return alleleCount;
}

inline float Bgen::getFrequency()
{
	return frequency;
}

inline float Bgen::getLLD()
{
	return LLD;
}


#endif
