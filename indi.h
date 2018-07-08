#ifndef INDI_H_
#define INDI_H_


#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <string>


class Indi
{
public: 
	Indi(const char* IdFileName, const char* positionFileName, int partNum, int totalNum);
	int getN();
	std::vector<int>* getIds();
	std::vector<int> BgenPositions;

private:
	// all the IDs that could be used as read in from a file
	std::vector<int> allIds;
	std::vector<int> allBgenPositions;

	// subset of IDs and corresponding positions in the BGEN file that are to be used here
	std::vector<int> Ids;
	int n; // number of individuals used in the simulation
};


inline std::vector<int>* Indi::getIds()
{
	return &Ids; 
}

inline int Indi::getN()
{
	return n;
}






#endif
