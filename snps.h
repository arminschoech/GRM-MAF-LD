#ifndef SNPS_H_
#define SNPS_H_


#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <string>


class Snps
{
public: 
	Snps(const char* fileName, int inIndiN, int inMinimalAlleleCount, bool lldInfoYN);
	bool useSnp(std::string rsId, std::string snpId);
	float getSnpLD(std::string rsId, std::string snpId);

private:
	int indiN; //number of individuals used in the frequency file
	int snpN; //number of snps in that file
	int minimalAlleleCount;
	std::vector<std::string> snpIdVec; 
	std::vector<int> alleleCountVec;
	std::vector<float> LDvec;

	int getAlleleCount(std::string rsId, std::string snpId);
	int binarySearch(std::vector<std::string>* vec, std::string word);
};








#endif
