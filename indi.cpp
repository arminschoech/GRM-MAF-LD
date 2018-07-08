#include "indi.h"



Indi::Indi(const char* IdFileName, const char* positionFileName, int partNum, int totalNum)
{
  // read in ID file
  std::ifstream inIdFile(IdFileName);
	if (!inIdFile.is_open()) {
    std::cerr << "ERROR: ID file could not be opened." << std::endl;
    exit(1);
  }

  std::string line;
  while (getline(inIdFile,line)) {
    allIds.push_back(atoi(line.c_str()));
  }
    
  int allIdsNum = allIds.size();
  inIdFile.close();


    // read in BGEN position file
	std::ifstream inPosFile(positionFileName);
	if (!inPosFile.is_open()) {
    std::cerr << "ERROR: BGEN position file could not be opened." << std::endl;
    exit(1);
  }

  while (getline(inPosFile,line)) {
    allBgenPositions.push_back(atoi(line.c_str()));
  }

  int allBgenPositionsNum = allBgenPositions.size();
  inPosFile.close();

  if (allBgenPositionsNum != allIdsNum) {
    std::cerr << "ERROR: BGEN position file and ID file do not have the same number of individuals." << std::endl;
    exit(1);
  }

  // create ID and position files for the individuals used in this simulation
  // taking all the individuals in the used ID file and taking the ith part out of n equally large parts
  // (just the last parts might have a couple more individuals if the number is not a multiple)
  int sizeOfMostPatches = allIdsNum / totalNum;
  int sizeOfLastPatch = allIdsNum - (totalNum-1) * sizeOfMostPatches;

  int sizeOfPatch;
  if(partNum == totalNum) {
    sizeOfPatch = sizeOfLastPatch;
  } else {
    sizeOfPatch = sizeOfMostPatches; 
  }

  int startNum = (partNum-1)*sizeOfMostPatches;
  for(int i = startNum; i < startNum+sizeOfPatch; i++) {
    Ids.push_back(allIds[i]);
    BgenPositions.push_back(allBgenPositions[i]);
  }

  n = sizeOfPatch;
}


