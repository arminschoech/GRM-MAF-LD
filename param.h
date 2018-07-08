#ifndef param_H_
#define param_H_


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
using std::cout;
using std::endl;
using std::cerr;

class Param 
{
public: 
	Param(int argc, char* argv[]);
	void printParameterSummary(int argc, char* argv[]);
	int getSeed();
	int getIndiPartNum();
	int getTotalIndiPartNum();
	int getSnpPartNum();
	int getTotalSnpPartNum();
	float getAlpha();
	float getTau();
	float getCausalFrac();
	float getHeritability();
	const char* getBgenFileName();
	const char* getFreqFileName();
	int getRefPanelNum();
	bool getSimPhenoYN();
	bool getSampledGenoYN();
	const char* getOutFileName();
	int getMBuffer();
	int getMinimalAlleleCount();
	int getMinimalAlleleCountRefPanel();
	const char* getIndividualIdFileName();
	const char* getBgenPositionFileName();
	bool getLldInfoYN();

private:
	int seed; 
	int indiPartNum;
       	int totalIndiPartNum; 
	int snpPartNum; 
	int totalSnpPartNum; 
	float alpha; 
	float tau; 
	float causalFrac; 
	float heritability;
	const char* bgenFileName; 
	const char* freqFileName;	
	int refPanelNum;
	bool simPhenoYN;
	bool sampledGenoYN;
	const char* outFileName;
	int mBuffer;
	int minimalAlleleCount;
	int minimalAlleleCountRefPanel;
	const char* individualIdFileName;
	const char* bgenPositionFileName;
	bool lldInfoYN;
	void setLldInfoYN();
};


inline int Param::getSeed()
{
	return seed;
}

inline int Param::getIndiPartNum()
{
	return indiPartNum;
}

inline int Param::getTotalIndiPartNum()
{
	return totalIndiPartNum;
}

inline int Param::getSnpPartNum()
{
	return snpPartNum;
}

inline int Param::getTotalSnpPartNum()
{
	return totalSnpPartNum;
}

inline float Param::getAlpha()
{
	return alpha;
}

inline float Param::getTau()
{
	return tau;
}

inline float Param::getCausalFrac()
{
	return causalFrac;
}

inline float Param::getHeritability()
{
	return heritability;
}

inline const char* Param::getBgenFileName()
{
	return bgenFileName;
}

inline const char* Param::getFreqFileName()
{
	return freqFileName;
}

inline int Param::getRefPanelNum()
{
	return refPanelNum;
}

inline bool Param::getSimPhenoYN()
{
	return simPhenoYN;
}

inline const char* Param::getOutFileName()
{
	return outFileName;
}

inline int Param::getMBuffer()
{
	return mBuffer;
}

inline int Param::getMinimalAlleleCount()
{
	return minimalAlleleCount;
}

inline int Param::getMinimalAlleleCountRefPanel()
{
	return minimalAlleleCountRefPanel;
}

inline bool Param::getSampledGenoYN()
{
	return sampledGenoYN;
}

inline const char* Param::getIndividualIdFileName()
{
	return individualIdFileName;
}

inline const char* Param::getBgenPositionFileName()
{
	return bgenPositionFileName;
}

inline bool Param::getLldInfoYN()
{
	return lldInfoYN;
}

#endif
