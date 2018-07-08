#include "bgen.h"

Bgen::Bgen(const char* fileName, int partNum, int totalNum, int inMinimalAlleleCount, Indi* inIndividuals, Snps* inSnpFile)
{
	bgenFile = fopen(fileName, "rb");
	if (bgenFile == NULL) {
    	std::cerr << "ERROR: BGEN file cannot be opened." << std::endl;
    	exit(1);
  	}

	individuals = inIndividuals;
	n = individuals->getN();

	snpFile = inSnpFile;

	minAlleleCount = inMinimalAlleleCount;
	loadSnpSuccess = true;


	// read in header information from BGEN file
  	unsigned int offset; fread(&offset, 4, 1, bgenFile); 
  	unsigned int H; fread(&H, 4, 1, bgenFile);
  	unsigned int uIntM; fread(&uIntM, 4, 1, bgenFile); mBgen = (int) uIntM;
  	unsigned int uIntN; fread(&uIntN, 4, 1, bgenFile); nBgen = (int) uIntN;
  	

  	// define start and end SNP of simulation
  	int batchM = mBgen / totalNum; // number of SNPs for all simulation unless it's the last part
  	startM = batchM * (partNum - 1);
  	if(partNum != totalNum) {
  		endM = batchM * partNum - 1;
  	} else {
  		endM = mBgen - 1; // last parts runs until very end (might have a couple more SNPs)
  	}
  	currM = 0;

  	// declare variables and buffer
  	fseek(bgenFile, offset+4, SEEK_SET);
  	zBuf = (unsigned char *) malloc(6*nBgen);
  	shortBuf = (unsigned short *) malloc(6*nBgen);

  	// initialize genotype vector
  	genVec = gsl_vector_float_alloc(n);

  	// skip SNPs until you reach startM
  	skipSnps(startM);
}


bool Bgen::useSnp()
{
	return (alleleCount >= minAlleleCount && alleleCount <= 2*n - minAlleleCount);
}


void Bgen::sampleNextNonZeroSnp()
{
	while(true) 
	{
		if (currM > endM) 
		{
    		loadSnpSuccess = false;
    		break;
    	}

		// read in genotype data
		int zLen = readSnpHeader();

		if(snpFile->useSnp(currRsId, currSnpId)) { //if UK10K allele frequency is high enough
			
			fread(zBuf, 1, zLen, bgenFile); //read in compressed data
			uLongf destLen = 6*nBgen;
    		if (uncompress((Bytef *) shortBuf, &destLen, zBuf, zLen) != Z_OK || destLen != 6 * (unsigned int) nBgen) {
      			std::cerr << "ERROR: uncompress() failed" << std::endl;
      			exit(1);
    		}

			sampleSnp();
    		currM++;

    		if(useSnp()) { //if UK Biobank frequency is high enough
    			LLD = snpFile->getSnpLD(currRsId, currSnpId);
    			break;
    		}

    	} else {
    		fseek(bgenFile, zLen, SEEK_CUR);
			currM++;
    	}
    }
}


void Bgen::mapNextNonZeroSnp()
{
	while(true) 
	{
		if (currM > endM) 
		{
    		loadSnpSuccess = false;
    		break;
    	}

		// read in genotype data
		int zLen = readSnpHeader();

		if(snpFile->useSnp(currRsId, currSnpId)) { //if UK10K allele frequency is high enough
			
			fread(zBuf, 1, zLen, bgenFile); //read in compressed data
			uLongf destLen = 6*nBgen;
    		if (uncompress((Bytef *) shortBuf, &destLen, zBuf, zLen) != Z_OK || destLen != 6 * (unsigned int) nBgen) {
      			std::cerr << "ERROR: uncompress() failed" << std::endl;
      			exit(1);
    		}

			mapSnp();
    		currM++;

    		if(useSnp()) { //if UK Biobank frequency is high enough
    			LLD = snpFile->getSnpLD(currRsId, currSnpId);
    			break;
    		}

    	} else {
    		fseek(bgenFile, zLen, SEEK_CUR);
			currM++;
    	}
    }
}


void Bgen::dosageNextNonZeroSnp()
{
	while(true) 
	{
		if (currM > endM) 
		{
    		loadSnpSuccess = false;
    		break;
    	}

		// read in genotype data
		int zLen = readSnpHeader();

		if(snpFile->useSnp(currRsId, currSnpId)) { //if UK10K allele frequency is high enough
			
			fread(zBuf, 1, zLen, bgenFile); //read in compressed data
			uLongf destLen = 6*nBgen;
    		if (uncompress((Bytef *) shortBuf, &destLen, zBuf, zLen) != Z_OK || destLen != 6 * (unsigned int) nBgen) {
      			std::cerr << "ERROR: uncompress() failed" << std::endl;
      			exit(1);
    		}

			dosageSnp();
    		currM++;

    		if(useSnp()) { //if UK Biobank frequency is high enough
    			LLD = snpFile->getSnpLD(currRsId, currSnpId);
    			break;
    		}

    	} else {
    		fseek(bgenFile, zLen, SEEK_CUR);
			currM++;
    	}
    }
}


int Bgen::readSnpHeader()
{
	unsigned int Nrow; fread(&Nrow, 4, 1, bgenFile);
    if (Nrow != (unsigned int) nBgen) {
      std::cerr << "ERROR: Nrow = " << Nrow << " does not match N = " << nBgen << std::endl;
      exit(1);
    }
    unsigned short LS; fread(&LS, 2, 1, bgenFile); 
    char snpID[1000]; fread(snpID, 1, LS, bgenFile); snpID[LS] = '\0'; currSnpId = std::string(snpID);
    unsigned short LR; fread(&LR, 2, 1, bgenFile); 
    char rsID[1000]; fread(rsID, 1, LR, bgenFile); rsID[LR] = '\0'; currRsId = std::string(rsID);
    unsigned short LC; fread(&LC, 2, 1, bgenFile); 
    char chrStr[1000]; fread(chrStr, 1, LC, bgenFile); chrStr[LC] = '\0';
    int chrom;
    if (sscanf(chrStr, "%d", &chrom)!=1 || !(1<=chrom&&chrom<=22)) {
      std::cerr << "ERROR: Invalid chrom (expecting integer 1-22): " << std::string(chrStr) << std::endl;
      exit(1);
    }
    unsigned int physpos; fread(&physpos, 4, 1, bgenFile); // cout << "physpos: " << physpos << endl;
    unsigned int LA; fread(&LA, 4, 1, bgenFile); // cout << "LA: " << LA << endl;
    char allele1[1000]; fread(allele1, 1, LA, bgenFile); allele1[LA] = '\0';
    unsigned int LB; fread(&LB, 4, 1, bgenFile); // cout << "LB: " << LB << endl;
    char allele0[1000]; fread(allele0, 1, LB, bgenFile); allele0[LB] = '\0';
    unsigned int zLen; fread(&zLen, 4, 1, bgenFile); // cout << "zLen: " << zLen << endl;

    return (int) zLen;
}


void Bgen::skipSnps(int skipM)
{
	for(int i = 0; i < skipM; i++)
	{
		// read in genotype data
		int zLen = readSnpHeader(); 
		fseek(bgenFile, zLen, SEEK_CUR);
		currM++;
	}
}


void Bgen::mapSnp()
{
	alleleCount = 0;
	unsigned short AA;
	unsigned short AB;
	unsigned short BB;

	for (int i = 0; i < n; i++) { 
		// look up the BGEN file position of the ith used individual
		// and get the genotype probabilities there
		AA = shortBuf[ 3*individuals->BgenPositions[i] ];
		AB = shortBuf[ 3*individuals->BgenPositions[i] + 1];
		BB = shortBuf[ 3*individuals->BgenPositions[i] + 2];

		if(AA > AB && AA > BB)
		{
			gsl_vector_float_set(genVec, i, 0.0);
		} 
		else if(AB > BB)
		{
			gsl_vector_float_set(genVec, i, 1.0); alleleCount += 1;
		}
		else
		{
			gsl_vector_float_set(genVec, i, 2.0); alleleCount += 2;
		}
    }
    frequency = (float) alleleCount / 2.0 / (float) n;
}


void Bgen::sampleSnp()
{
	int intAA;
	int intAB;
	int intBB;
	int genoSum;

	int randomNumber;

	alleleCount = 0;

	for(int i = 0; i < n; i++)
	{
		// look up the BGEN file position of the ith used individual
		// and get the genotype probabilities there
		intAA = (int) shortBuf[ 3*individuals->BgenPositions[i] ];
		intAB = (int) shortBuf[ 3*individuals->BgenPositions[i] + 1];
		intBB = (int) shortBuf[ 3*individuals->BgenPositions[i] + 2];
		genoSum = intAA + intAB + intBB;

		if(genoSum == 0) // zero probability for all three genotypes throws error 
		{
			std::cerr << "ERROR: Individual " << i << " at SNP " << currM << " has 0 probability for all three genotypes!" << std::endl; 
			exit(1);
		}

		if(genoSum == intAA) // if prob(AA) = 1
		{
			gsl_vector_float_set(genVec, i, 0.0); 
			continue;
		}

		if(genoSum == intAB) // if prob(AB) = 1
		{
			gsl_vector_float_set(genVec, i, 1.0); alleleCount += 1; 
			continue;
		}

		if(genoSum == intBB) // if prob(BB) = 1
		{
			gsl_vector_float_set(genVec, i, 2.0); alleleCount += 2; 
			continue;
		}

		// if there is indeed a chance for two different genotypes
		randomNumber = gsl_rng_uniform_int(ran, genoSum);
		
		if(randomNumber < intAA)
		{
			gsl_vector_float_set(genVec, i, 0.0); 

		} 
		else if(randomNumber >= intAA + intAB)
		{
			gsl_vector_float_set(genVec, i, 2.0); alleleCount += 2;
		}
		else
		{
			gsl_vector_float_set(genVec, i, 1.0); alleleCount += 1;
		}
	}
	frequency = (float) alleleCount / 2.0 / (float) n;
}


void Bgen::dosageSnp()
{
	float AA;
	float AB;
	float BB;
	float genoSum;

	double dosageSum = 0;
	float dosage; 

	for(int i = 0; i < n; i++)
	{
		// look up the BGEN file position of the ith used individual
		// and get the genotype probabilities there
		AA = float(shortBuf[ 3*individuals->BgenPositions[i] ]);
		AB = float(shortBuf[ 3*individuals->BgenPositions[i] + 1]);
		BB = float(shortBuf[ 3*individuals->BgenPositions[i] + 2]);
		genoSum = AA + AB + BB;
		
		dosage = (AB + 2.0*BB) / genoSum;
		gsl_vector_float_set(genVec, i, dosage); 
		dosageSum += dosage; 
	}
	alleleCount = int(dosageSum + 0.5);
	frequency = dosageSum / 2.0 / (float) n;
}



void Bgen::printGenProb(int maxN)
{
	std::cout << "Bgen probability matrix:" << std::endl;
	for(int i = 0; i < maxN; i++)
	{
		std::cout << shortBuf[3*i]/32768.0 << '\t' << shortBuf[3*i+1]/32768.0 << '\t' << shortBuf[3*i+2]/32768.0 << std::endl;
	}
}

void Bgen::printGenProb()
{
	std::cout << "Bgen probability matrix:" << std::endl;
	for(int i = 0; i < n; i++)
	{
		std::cout << shortBuf[3*i]/32768.0 << '\t' << shortBuf[3*i+1]/32768.0 << '\t' << shortBuf[3*i+2]/32768.0 << std::endl;
	}
}

void Bgen::printGenVec(int maxN)
{
	std::cout << "Genotype vector:" << std::endl;
	for(int i = 0; i < maxN; i++) {
		std::cout << gsl_vector_float_get(genVec, i) << std::endl;
	}
	std::cout << std::endl;
}












