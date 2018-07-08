#include "geno.h"
#include "grm.h"
#include "pheno.h"
#include "bgen.h"
#include "indi.h"
#include "snps.h"
#include "random.h"
#include "param.h"

#include "omp.h" 
using std::cout;
using std::endl;

// declare global random number generator object
gsl_rng* ran;


int main(int argc, char* argv[])
{
	// making sure that only one core is used
	omp_set_num_threads(1);
	cout << endl << "+----------------------------------------------------------+";
	cout << endl << "| Tool for fast frequency and LD dependent GRM calculation |";
	cout << endl << "+----------------------------------------------------------+" << endl << endl;

	// read in parameters and print parameter summary 
	Param par(argc, argv);
	par.printParameterSummary(argc, argv);

	// seed random number generator:
	ran = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(ran, par.getSeed());

	// initialize objects
	Indi individuals(par.getIndividualIdFileName(), par.getBgenPositionFileName(), par.getIndiPartNum(), par.getTotalIndiPartNum());
	//TODO: potentially add default for ID and BGENposition file
	Snps snpFile(par.getFreqFileName(), par.getRefPanelNum(), par.getMinimalAlleleCountRefPanel(), par.getLldInfoYN());
	Bgen bgenLink(par.getBgenFileName(), par.getSnpPartNum(), par.getTotalSnpPartNum(), par.getMinimalAlleleCount(), &individuals, &snpFile);
	Geno genotypes(individuals.getN(), par.getMBuffer(), par.getAlpha(), par.getTau());
	Pheno phenotypes(individuals.getN());
	int grmN = individuals.getN();
	if(par.getSimPhenoYN()) grmN = 10; //avoid large memory cost when GRM is not needed
	Grm myGrm(grmN, par.getAlpha()); 

	// reading in genotypes and calculating GRM/phenotypes
	cout << endl << "Processing genotype data ..." << endl; 
	while(bgenLink.loadSnpSuccess) // as long as the last buffer genotype matrix was fully filled
	{
		genotypes.resetGenotypeMatrix();

		for(int j = 0; j < par.getMBuffer(); j++)
		{
			if(par.getSampledGenoYN()) 
				bgenLink.sampleNextNonZeroSnp(); //randomly drawn from BGEN prob.
			else
				bgenLink.dosageNextNonZeroSnp(); //using genotype dosages

			if(bgenLink.loadSnpSuccess) { // if last SNP was successfully loaded
				genotypes.readInGenVec(&bgenLink);
			} else { // if there are no more SNPs to read in
				genotypes.fillRest(); // fill rest of matrix with zeros
				break;
			}
		}
		if(par.getSimPhenoYN())
			phenotypes.addGenEffects(&genotypes, par.getCausalFrac());
		else
			myGrm.addGenotypes(&genotypes);
	}
	if(!par.getSimPhenoYN()) myGrm.divideBySnpNumber();

	if(par.getSimPhenoYN()) {
		phenotypes.writePhenoFile(par.getHeritability(), par.getOutFileName(), &individuals);
		cout << endl << "Writing phenotype values to '" << par.getOutFileName() << ".phen'." << endl;
	} else {
		cout << endl << "Writing binary GRM files " << par.getOutFileName() << ".grm.bin, " << par.getOutFileName() << ".grm.N.bin and " << par.getOutFileName() << ".grm.id." << endl << endl;
		myGrm.writeGrmFiles(par.getOutFileName(), &individuals);
	}
	cout << "Analysis finished." << endl;
	return 0;
}


