#include "param.h"

Param::Param(int argc, char* argv[])
{
	// default parameters
	seed = 0; 
	indiPartNum = 1;
	totalIndiPartNum = 1;
	snpPartNum = 1;
	totalSnpPartNum = 1;
	tau = 0.0;
	causalFrac = 1.0;
	heritability = 1.0; 
	simPhenoYN = false; 
	sampledGenoYN = false;
	outFileName = "out";
	refPanelNum = 3567;
	lldInfoYN = false;

	//fixed parameters
	mBuffer = 1000;
	minimalAlleleCount = 2;
	minimalAlleleCountRefPanel = 5;

	//checks that parameters have been assigned
	bool seedBool = false;
	bool alphaBool = false;
	bool bgenBool = false;
	bool freqBool = false;
	bool idBool = false;
	bool bgenPosBool = false;
	bool h2bool = false;

	for(int i = 1; i < argc; i++)
	{
		if (!strncmp(argv[i], "--seed", 100)) {
			seed = atoi(argv[++i]);
			seedBool = true;
		}
		else if (!strncmp(argv[i], "--alpha", 100)) {
			alpha = atof(argv[++i]);
			alphaBool = true;
		}
		else if (!strncmp(argv[i], "--tau", 100)) 
			tau = atof(argv[++i]);
		else if (!strncmp(argv[i], "--fraction_causal", 100)) 
			causalFrac = atof(argv[++i]);
		else if (!strncmp(argv[i], "--h2", 100)) {
			heritability = atof(argv[++i]);
			h2bool = true;
		}
		else if (!strncmp(argv[i], "--bgen_file", 100)) {
			bgenFileName = argv[++i];
			bgenBool = true;
		}
		else if (!strncmp(argv[i], "--freq_file", 100)) { 
			freqFileName = argv[++i];
			freqBool = true;
		}
		else if (!strncmp(argv[i], "--ref_panel_num", 100))  
			refPanelNum = atoi(argv[++i]);
		else if (!strncmp(argv[i], "--out", 100)) 
			outFileName = argv[++i];
		else if (!strncmp(argv[i], "--individual_partition", 100)) {
			indiPartNum = atoi(argv[++i]);
			totalIndiPartNum = atoi(argv[++i]);
		}
		else if (!strncmp(argv[i], "--snp_partition", 100)) {
			snpPartNum = atoi(argv[++i]);
			totalSnpPartNum = atoi(argv[++i]);
		}
		else if (!strncmp(argv[i], "--sim_pheno", 100)) 
			simPhenoYN = true;
		else if (!strncmp(argv[i], "--sampled_genotypes", 100)) 
			sampledGenoYN = true;
		else if (!strncmp(argv[i], "--id_file", 100)) { 
			individualIdFileName = argv[++i];	
			idBool = true;
		}
		else if (!strncmp(argv[i], "--bgen_pos_file", 100)) {
			bgenPositionFileName = argv[++i];
			bgenPosBool = true;
		}
		else {
			cerr << "ERROR: '" << argv[i] << "' is not a valid input flag." << endl;
			exit(1);
		}
	}

	// raise errors when necessary parameters are not specified
	if(!alphaBool) {
		cerr << "ERROR: '--alpha' flag has to be specified." << endl;
		exit(1);
	}
	if(!bgenBool) {
		cerr << "ERROR: '--bgen_file' flag has to be specified." << endl;
		exit(1);
	}
	if(!freqBool) {
		cerr << "ERROR: '--freq_file' flag has to be specified." << endl;
		exit(1);
	}
	if(simPhenoYN && !seedBool) {
		cerr << "ERROR: when simulating phenotypes, the '--seed' flag has to be specified." << endl;
		exit(1);
	} else if(sampledGenoYN && !seedBool) {
		cerr << "ERROR: when using sampled genotypes, the '--seed' flag has to be specified." << endl;
		exit(1);
	}
	if(simPhenoYN && !h2bool) {
		cerr << "ERROR: when simulating phenotypes, the '--h2' flag has to be specified." << endl;
		exit(1);
	}
	if(!bgenPosBool) {
		cerr << "ERROR: '--bgen_pos_file' flag has to be specified." << endl;
		exit(1);
	}
	if(!idBool) {
		cerr << "ERROR: '--id_file' flag has to be specified." << endl;
		exit(1);
	}

	setLldInfoYN();
	if(tau != 0.0 && !lldInfoYN) {
		std::cerr << "ERROR: If the tau parameter is set to a value other than zero, LLD values have to be added as an additional column in the frequency file." << std::endl;
		exit(1);
	}
}


void Param::printParameterSummary(int argc, char* argv[]) {
	cout << endl << "Options:";
	for(int i = 1; i < argc ; i++) {
		if(argv[i][0] == '-' && argv[i][1] == '-') 
			cout << endl << "  ";	
		cout << argv[i] << ' ';
	}
	cout << endl;
}

void Param::setLldInfoYN() {
	std::ifstream frqFile(freqFileName);
	if (!frqFile.is_open()) {
		std::cerr << "ERROR: UK10K frequency file could not be opened." << std::endl;
		exit(1);
	}

	std::string word;
	frqFile >> word; // chromosome
	frqFile >> word; // snp ID
	frqFile >> word; // primary allele
	frqFile >> word; // secondary allele
	frqFile >> word; // freq
	frqFile >> word; // NCHROBS
	frqFile >> word; // LLD??

	if(word == "LLD") 
		lldInfoYN = true;
}





