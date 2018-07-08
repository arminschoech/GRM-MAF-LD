#ifndef PHENO_H_
#define PHENO_H_


#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <gsl/gsl_vector.h>

#include "geno.h"
#include "random.h"




class Pheno // stores and manipulates a vector of phenotype values for all individuals
{
public:
	Pheno(int inN);
	void addGenEffects(Geno*, float causalFrac);
	void writePhenoFile(double h2, std::string fileName, Indi* indiFile);

	void printPhenotypes(); // print all phenotypes to terminal
	void printPhenotypes(int maxN); // print the first Nmax phenotypes to terminal

private:
	int n;
	gsl_vector_float* pPhenotypes; // only includes genetic effects but not environmental effects

	gsl_vector_float* simulateEffectVector(int snpN, float causalFrac);

	double getGeneticVariance();

};












#endif
