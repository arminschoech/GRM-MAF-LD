#include "snps.h"


Snps::Snps(const char* fileName, int inIndiN, int inMinimalAlleleCount, bool lldInfoYN)
{
	indiN = inIndiN;
	minimalAlleleCount = inMinimalAlleleCount;

	std::ifstream frqFile(fileName);
	if (!frqFile.is_open()) {
    	std::cerr << "ERROR: UK10K frequency file could not be opened." << std::endl;
    	exit(1);
 	}
 	
 	// get number of lines
 	snpN = -1; // -1 to account for header line
	std::string line;
	while (std::getline(frqFile , line))
		snpN++;
	
	snpIdVec.resize(snpN);
	alleleCountVec.resize(snpN);
	LDvec.resize(snpN);	
	for(int i = 0; i < snpN; i++) {
		LDvec[i] = 0;
	}

	// put stream pointer back to beginning
	frqFile.clear();
	frqFile.seekg(0, std::ios_base::beg);

	// define 
 	std::string stringFreq;
 	std::string stringLD;
 	std::string trash;
 	float maf; // minor allele frequency

 	// reading in header
 	frqFile >> trash; // chromosome 	 		
 	frqFile >> trash; // snp ID
 	frqFile >> trash; // primary allele
 	frqFile >> trash; // secondary allele
 	frqFile >> trash; // freq
 	frqFile >> trash; // NCHROBS
	if(lldInfoYN) {
 		frqFile >> trash; // LLD
	}

 	for(int i = 0; i < snpN; i++) {
 		frqFile >> trash; // chromosome 	 		
 		frqFile >> snpIdVec[i];
 		frqFile >> trash; // primary allele
 		frqFile >> trash; // secondary allele
 		frqFile >> stringFreq; maf = atof(stringFreq.c_str()); alleleCountVec[i] = int(maf * 2.0 * float(indiN) + 0.5);
 		frqFile >> trash; // NCHROBS
		if(lldInfoYN) {
 			frqFile >> stringLD; LDvec[i] = atof(stringLD.c_str());
		}
 	}

	frqFile.close();
}


bool Snps::useSnp(std::string rsId, std::string snpId)
{
	int count = getAlleleCount(rsId, snpId);
	return (count >= minimalAlleleCount);
}

int Snps::getAlleleCount(std::string rsId, std::string snpId)
{
	int matchPosition = binarySearch(&snpIdVec, rsId);
	if(matchPosition == -1) matchPosition = binarySearch(&snpIdVec, snpId);

	if(matchPosition == -1) {
		return -1;
	} else {
		return alleleCountVec[matchPosition];
	}
}


float Snps::getSnpLD(std::string rsId, std::string snpId)
{
	int matchPosition = binarySearch(&snpIdVec, rsId);
	if(matchPosition == -1) matchPosition = binarySearch(&snpIdVec, snpId); // since we have to try both possible identifiers

	if(matchPosition == -1) {
		return 99.0;
	} else {
		return LDvec[matchPosition];
	}
}


int Snps::binarySearch(std::vector<std::string>* vec, std::string word)
{
	int matchPosition = -1;
	int lower = -1;
	int upper = vec->size();
	int check;
	std::string currWord;

	while(upper - lower >= 2) {
		check = (upper - lower)/2 + lower;
		currWord = vec->at(check);
		if(word > currWord) {
			lower = check;
		} else if(word < currWord) {
			upper = check;
		} else {
			matchPosition = check;
			break;
		}
	}

	return matchPosition;
}



