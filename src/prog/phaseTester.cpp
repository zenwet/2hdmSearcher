#include "EDM.h"
#include "HelpFunctions.h"
#include "SM.h"
#include "THDM.h"

#include <complex>
#include <iostream>
#include <cmath>

#include <time.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace THDME;
using namespace std;

void phaseTester(double betaAngle, double mOneTwo, double lambdaOne, double lambdaTwo, double lambdaThree, double lambdaFour, double lambdaFive, Eigen::Matrix3cd rho_lambdas){
	
	time_t rawtime;
	
	time (&rawtime);
	
	std::string curTime = ctime(&rawtime);
	
	std::string currTime = curTime.substr(0, curTime.size() - 1);
	
	currTime.erase(0,4);
	
	std::cout << "phaseTester initialised on " << currTime << "\n\n";
	
	FileSystem _fileManager;
	std::string _outDir = "output/data/";
	_fileManager.set_directory(_outDir);
	
	std::string _csvOut = currTime + "_phaseTester" + ".csv";
	_fileManager.create_file(_csvOut);
	
	std::vector<std::string> _header;
	
	_header.push_back("Lambda5 phase");
	
	for(int h = 0; h < 3; ++h){
		std::string currMh = "mh_" + std::to_string(h + 1);
		_header.push_back(currMh);
	}
	
	_header.push_back("mh_C");
	
	_header.push_back("eEDM");
	_header.push_back("muAMM");
	
	_header.push_back("1-loop eEDM contr");
	_header.push_back("1-loop muAMM contr");
	
	_header.push_back("tan(beta)");
	_header.push_back("M12");
	for(int lamb = 1; lamb <=4; ++lamb){
		std::string curren = std::to_string(lamb);
		std::string currLambda = "Lambda" + curren;
		_header.push_back(currLambda);
	}
	
	for(int lamb = 5; lamb <=7; ++lamb){
		std::string curren = std::to_string(lamb);
		std::string currLambda = "Lambda" + curren;
		_header.push_back("Re(" + currLambda + ")");
		_header.push_back("Im(" + currLambda + ")");
	}
	
	_header.push_back("rhoL mag factor");
	
	for(int row = 0; row < 3; ++row){
		for(int colum = 0; colum < 3; ++ colum){
		std::string currRho;
		currRho = "lambdaL_(" + std::to_string(row + 1) + "," + std::to_string(colum + 1) + ")";
		_header.push_back(currRho);
		}
	}

	_fileManager.add_line(_csvOut, _header);
	
	_header.clear();
	
	std::complex<double> _i(0., 1.);

	int resolution = 200;
	
	
	for(int tems = 0; tems <= resolution; ++tems){
	
		SM _sm;
	
		THDM _thdm(_sm);
	
		Base_generic gen;
		
		std::vector<double> outputLine;
		
		double _currPhase = M_PI * tems / (1. * resolution);
		std::complex<double> _currPhaseComp = std::exp(_i * _currPhase);
		
		gen.beta = betaAngle;
		gen.M12 = mOneTwo;
		gen.Lambda1 = lambdaOne;
		gen.Lambda2 = lambdaTwo;
		gen.Lambda3 = lambdaThree;
		gen.Lambda4 = lambdaFour;
		gen.Lambda5 = lambdaFive * _currPhaseComp;
		gen.Lambda6 = std::complex<double>(0., 0.);
		gen.Lambda7 = std::complex<double>(0., 0.);
	
		outputLine.push_back(std::real(_currPhase));
		
		Eigen::Matrix3cd rhoU;
		Eigen::Matrix3cd rhoD;
		Eigen::Matrix3cd rhoL;
		
		double _v2 = _thdm.get_v2();
		
		FermionSector up = UP;
		FermionSector down = DOWN;
		FermionSector lepton = LEPTON;
		
		for(int rowpos = 0; rowpos < 3; ++rowpos){
			for(int colpos = 0; colpos < 3; ++colpos){
				if(rowpos == colpos){
					rhoU(rowpos, colpos) = 	std::complex<double>(_thdm.get_mf(up, rowpos + 1) * std::sqrt(2. / _v2), 0.);
					rhoD(rowpos, colpos) = std::complex<double>(_thdm.get_mf(down, rowpos + 1) * std::sqrt(2. / _v2), 0.);
					rhoL(rowpos, colpos) = std::complex<double>(_thdm.get_mf(lepton, rowpos + 1) * std::sqrt(2. / _v2), 0.);
				} else{
					rhoL(rowpos, colpos) = std::complex<double>(rho_lambdas(rowpos, colpos) * std::sqrt(_thdm.get_mf(lepton, rowpos + 1) * _thdm.get_mf(lepton, colpos + 1) * 2. / _v2));
				}
			}
		}
		
		_thdm.set_param_gen(gen);
		
		_thdm.set_yukawa_manual(rhoU, rhoD, rhoL);
		
		vector<double> _higgsMasses = _thdm.get_higgs_treeLvl_masses();
		
		for(int higgs = 0; higgs <= 3; ++higgs){
			outputLine.push_back(std::real(_higgsMasses[higgs]));
		}
		EDM _edm(_thdm);
		
		double _eEDM = _edm.get_electron_edm();
		double _aAMM = _edm.get_muon_amm();

		double _eEDM_1loop = 0.;
		double _aAMM_1loop = 0.;
				
		std::vector<double> _contribsEDM = _edm.get_de();
		std::vector<double> _contribsAMM = _edm.get_amu();
				
		for(int h = 0; h < 3; ++h){
			_eEDM_1loop += _contribsEDM[1 + 12*h];
			_aAMM_1loop += _contribsAMM[1 + 12*h];
		}

		outputLine.push_back(std::real(_eEDM));
		outputLine.push_back(std::real(_aAMM));
		
		outputLine.push_back(std::real(_eEDM_1loop));
		outputLine.push_back(std::real(_aAMM_1loop));

		if(tems == 0){
			outputLine.push_back(std::real(std::tan(gen.beta)));
			outputLine.push_back(std::real(gen.M12));
			outputLine.push_back(std::real(gen.Lambda1));
			outputLine.push_back(std::real(gen.Lambda2));
			outputLine.push_back(std::real(gen.Lambda3));
			outputLine.push_back(std::real(gen.Lambda4));
			outputLine.push_back(std::real(gen.Lambda5));
			outputLine.push_back(std::imag(gen.Lambda5));
			outputLine.push_back(std::real(gen.Lambda6));
			outputLine.push_back(std::imag(gen.Lambda6));
			outputLine.push_back(std::real(gen.Lambda7));
			outputLine.push_back(std::imag(gen.Lambda7));
			
			for(int rowpos = 0; rowpos < 3; ++rowpos){
				for(int colpos = 0; colpos < 3; ++colpos){
					outputLine.push_back(std::real(rho_lambdas(rowpos, colpos)));
				}
			}
		}
		_fileManager.add_line(_csvOut, outputLine);

		_higgsMasses.clear();
		
		outputLine.clear();
		
		if( (tems % 10) == 0){
			std::cout << "Finished " << tems << " out of " << resolution << " point studies.\n";
		}
		
		
	}
	
	std::cout << "\n\nPhase study completed.\n\n";
}


int main(){
	
	double lamb1 = 0.;
	double lamb2 = 0.;
	double lamb3 = 0.;
	double lamb4 = 0.;
	double lamb5 = 0.;
	double masscross = 0.;
	double tabeta = 0.;
	
	std::cout << "Welcome to the phaseTester program.\nThis program tests the variation of the Higgs boson masses and the contribution to d_e and a_mu for a 2HDM with a soft Z_2 symmetry-breaking potential.\n\n";
	
	std::cout << "The potential is written in the generic basis. \nPlease write your potential.\n";
	
	std::cout << "lambda_1 = ";
	std::cin >> lamb1;
	std::cout << "lambda_2 = ";
	std::cin >> lamb2;
	std::cout << "lambda_3 = ";
	std::cin >> lamb3;
	std::cout << "lambda_4 = ";
	std::cin >> lamb4;
	std::cout << "|lambda_5| = ";
	std::cin >> lamb5;
	std::cout << "|m_12^2| = ";
	std::cin >> masscross;
	std::cout << "tan(beta) = ";
	std::cin >> tabeta;

	Eigen::Matrix3cd rhoEll;

	rhoEll.setZero();


	std::cout << "\n\n" << "Next, you will set the rho^L matrix coefficients. Using the Cheng-Sher-parametrisation, you will set the individual (real) lambdas. Note that the other rho matrices are assumed to be identical to their kappa counterparts.\n\n";
	std::cout << "lambda_e,e = ";
	std::cin >> rhoEll(0,0);
	std::cout << "lambda_e,mu = ";
	std::cin >> rhoEll(0,1);
	std::cout << "lambda_e,tau = ";
	std::cin >> rhoEll(0,2);
	std::cout << "lambda_mu,e = ";
	std::cin >> rhoEll(1,0);
	std::cout << "lambda_mu,mu = ";
	std::cin >> rhoEll(1,1);
	std::cout << "lambda_mu,tau = ";
	std::cin >> rhoEll(1,2);
	std::cout << "lambda_tau,e = ";
	std::cin >> rhoEll(2,0);
	std::cout << "lambda_tau,mu = ";
	std::cin >> rhoEll(2,1);
	std::cout << "lambda_tau,tau = ";
	std::cin >> rhoEll(2,2);
	std::cout << "\n";
	
	phaseTester(std::atan(tabeta), masscross, lamb1, lamb2, lamb3, lamb4, lamb5, rhoEll);
	
	std::cout << "Program finished. \n\n";
}