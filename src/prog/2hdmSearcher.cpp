#include "EDM.h"
#include "HelpFunctions.h"
#include "SM.h"
#include "THDM.h"

#include <complex>
#include <iostream>
#include <random>
#include <chrono>
#include <cmath>

#include <time.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

std::complex<double> _i(0., 1.);

using namespace THDME;
using namespace std;


//sign function defined for the sake of randomly generating a sign that can be biased towards negative or positive. also necessary when varying the square of a value linearly, rather than the value
template <typename T> inline constexpr
int signum(T x, std::false_type is_signed) {
    return T(0) < x;
}
template <typename T> inline constexpr
int signum(T x, std::true_type is_signed) {
    return (T(0) < x) - (x < T(0));
}
template <typename T> inline constexpr
int signum(T x) {
    return signum(x, std::is_signed<T>());
}


//function for searching for 2HDMs of type I,II,X, or Y
void typeSearcher(std::string thdmtype, double mintanbeta, double maxtanbeta, double minmonetwo, double maxonetwo, double lambda1, double lambda2, double lambda3, double lambda4, 
	double lambda5, double lam5minphase, double lam5maxphase, double amu, double de, int datapoints){
	
	//defining time string for output file
	time_t rawtime;
	time (&rawtime);
	std::string curTime = ctime(&rawtime);
	std::string currTime = curTime.substr(0, curTime.size() - 1);
	currTime.erase(0,4);
	
	std::mt19937_64 rng;
    // initialize random number generator with time-dependent seed
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    std::uniform_real_distribution<double> unif(0, 1);
	std::uniform_real_distribution<double> sign(-1,1);
	std::uniform_real_distribution<double> beta(std::exp(std::atan(mintanbeta)), std::exp(std::atan(maxtanbeta))); //varies of exp(beta) linearly, effectively varying beta logarithmically
	std::uniform_real_distribution<double> m12(minmonetwo * minmonetwo, maxonetwo * maxonetwo); //generates square of |m^2_12| linearly, so as to vary neutral boson masses linearly 
	std::uniform_real_distribution<double> lam5phase(lam5minphase, lam5maxphase); //generates phase of lambda_5 linearly
    // ready to generate random numbers
	
	FileSystem _fileManager;
	std::string _outDir = "output/data/";
	_fileManager.set_directory(_outDir);
	
	std::string _csvCheck = currTime + "_2hdmTypeSuccessInit" + ".csv";
	_fileManager.create_file(_csvCheck);
	
	std::string _csvSubSuc = currTime + "_2hdmTypeEdmAmmSuccess" + ".csv";
	_fileManager.create_file(_csvSubSuc);
	
	std::vector<std::string> _header;
	_header.push_back("tan(beta)");
	_header.push_back("M12");
	for(int lamb = 1; lamb <=7; ++lamb){
		std::string curren = std::to_string(lamb);
		std::string currLambda = "Lambda" + curren;
		_header.push_back(currLambda);
	}
	
	_header.push_back("a^L");
	
	for(int row = 0; row < 3; ++row){
		for(int colum = 0; colum < 3; ++ colum){
		std::string currRho;
		currRho = "lambdaL_(" + std::to_string(row + 1) + "," + std::to_string(colum + 1) + ")";
		_header.push_back(currRho);
		}
	}
	
	for(int h = 0; h < 3; ++h){
		std::string currMh = "mh_" + std::to_string(h + 1);
		_header.push_back(currMh);
	}
	
	_header.push_back("mh_C");
	
	_header.push_back("eEDM");
	_header.push_back("aAMM");

	_header.push_back("1-loop EDM contribution");
	_header.push_back("1-loop AMM contribution");

	_fileManager.add_line(_csvCheck, _header);
	_fileManager.add_line(_csvSubSuc, _header);
	
	_header.clear();
	
	int totalSuccEDM = 0;
	
	int totalSubSuccDeAmm = 0;
	
	for(int tems = 1; tems <= datapoints; ++tems){
		
		SM _sm;
	
		THDM _thdm(_sm);
	
		Base_generic gen;
		
		std::vector<complex<double>> outputLine;
		
		gen.beta = std::log(beta(rng));
		double monetwo = std::sqrt(m12(rng));
		gen.M12 = std::complex<double>(monetwo, 0.);
		gen.Lambda1 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda1 * lambda1);
		gen.Lambda2 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda2 * lambda2);
		gen.Lambda3 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda3 * lambda3);
		gen.Lambda4 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda4 * lambda4);
		std::complex<double> _randPhase_5 = std::exp(_i * lam5phase(rng));
		gen.Lambda5 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda5 * lambda5) * _randPhase_5;
		gen.Lambda6 = std::complex<double>(0., 0.);
		gen.Lambda7 = std::complex<double>(0., 0.); 
	
		_thdm.set_param_gen(gen);
		
		double rhoLmagFac = 0;
		
		if(thdmtype == "I"){
			_thdm.set_yukawa_type(TYPE_I);
			rhoLmagFac = 1. / std::tan(gen.beta);
		} else if(thdmtype == "II"){
			_thdm.set_yukawa_type(TYPE_II);
			rhoLmagFac = -std::tan(gen.beta);
		} else if(thdmtype == "Y" or thdmtype == "III"){
			_thdm.set_yukawa_type(TYPE_III);
			rhoLmagFac = 1. / std::tan(gen.beta);
		} else if(thdmtype == "X" or thdmtype == "IV" or thdmtype == "IIII"){
			_thdm.set_yukawa_type(TYPE_IV);
			rhoLmagFac = -std::tan(gen.beta);
		}
		
		outputLine.push_back(std::tan(gen.beta));
		outputLine.push_back(gen.M12);
		outputLine.push_back(gen.Lambda1);
		outputLine.push_back(gen.Lambda2);
		outputLine.push_back(gen.Lambda3);
		outputLine.push_back(gen.Lambda4);
		outputLine.push_back(gen.Lambda5);
		outputLine.push_back(gen.Lambda6);
		outputLine.push_back(gen.Lambda7);
		
		outputLine.push_back(rhoLmagFac);
		
		Eigen::Matrix3cd rho_lambdas;
		
		rho_lambdas.setZero();
		
		for(int rowpos = 0; rowpos < 3; ++rowpos){
			for(int colpos = 0; colpos < 3; ++colpos){
				if(rowpos == colpos){
					rho_lambdas(rowpos, colpos) = 1.;
				}
				outputLine.push_back(rho_lambdas(rowpos, colpos));
			}
		}
		
		vector<double> _higgsMasses = _thdm.get_higgs_treeLvl_masses();
		
		for(int higgs = 0; higgs <= 3; ++higgs){
			outputLine.push_back(_higgsMasses[higgs]);
		}
		
		
		if( (120 < _higgsMasses[0]) && (_higgsMasses[0] < 130)){
			std::cout << "Higgs mass between 120 and 130 GeV\n";
			
			EDM _edm(_thdm);

			if(_edm.is_initialized()){
				
				std::cout << "EDM type object initialised\n";
				
				std::cout << tems << " out of " << datapoints << " random tests done\n";
				
				if( (_thdm.is_perturbative() ) && (_thdm.is_unitary() ) && (_thdm.is_stable() ) ){
				
					totalSuccEDM += 1;
				
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
				
					outputLine.push_back(_eEDM);
					outputLine.push_back(_aAMM);
				
					outputLine.push_back(_eEDM_1loop);
					outputLine.push_back(_aAMM_1loop);
					
					_fileManager.add_line(_csvCheck, outputLine);
				
					if( (_eEDM < de) && (amu < _aAMM)){
						_fileManager.add_line(_csvSubSuc, outputLine);
						std::cout << "d_e and a_mu fulfilled input conditions\n";
					
						totalSubSuccDeAmm += 1;
					}
					
				}
			}
		}
		
		_higgsMasses.clear();
		
		outputLine.clear();
		
	}
	
	std::cout << "\n\n\nFinalised\n\n" << "Successfully initialised " << totalSuccEDM << " EDMs, out of which " << totalSubSuccDeAmm << " were found in interesting regions for d_e and a_mu.\n\n\n";
}



void alignSearcher(double mintanbeta, double maxtanbeta, double minmonetwo, 
		double maxonetwo, double lambda1, double lambda2, double lambda3, double lambda4, double lambda5, double lam5minphase, double lam5maxphase, 
		double lambda6, double lam6minphase, double lam6maxphase, double lambda7, double lam7minphase, double lam7maxphase, 
		double upalmin, double upalmax, double upphasemin, double upphasemax, double downalmin, double downalmax, double downphasemin, double downphasemax, 
		double leptalmin, double leptalmax, double leptphasemin, double leptphasemax, double amu, double de, int datapoints){
	
	time_t rawtime;
	time (&rawtime);
	std::string curTime = ctime(&rawtime);
	std::string currTime = curTime.substr(0, curTime.size() - 1);
	currTime.erase(0,4);
	
	std::mt19937_64 rng;
    // initialize random number generator with time-dependent seed
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    std::uniform_real_distribution<double> unif(0, 1);
	std::uniform_real_distribution<double> sign(-1,1);
	std::uniform_real_distribution<double> beta(std::exp(std::atan(mintanbeta)), std::exp(std::atan(maxtanbeta)));
	std::uniform_real_distribution<double> m12(minmonetwo * minmonetwo, maxonetwo * maxonetwo);
	std::uniform_real_distribution<double> lam5phase(lam5minphase, lam5maxphase);
	std::uniform_real_distribution<double> lam6phase(lam5minphase, lam5maxphase);
	std::uniform_real_distribution<double> lam7phase(lam5minphase, lam5maxphase);
	std::uniform_real_distribution<double> upInter(upalmin, upalmax);
	std::uniform_real_distribution<double> upPhase(upphasemin, upphasemax);
	std::uniform_real_distribution<double> downInter(downalmin, downalmax);
	std::uniform_real_distribution<double> downPhase(downphasemin, downphasemax);
	std::uniform_real_distribution<double> leptonInter(leptalmin, leptalmax);
	std::uniform_real_distribution<double> leptonPhase(leptphasemin, leptphasemax);
    // ready to generate random numbers
	
	FileSystem _fileManager;
	std::string _outDir = "output/data/";
	_fileManager.set_directory(_outDir);
	
	std::string _csvCheck = currTime + "_2hdmAlignedSuccessInit" + ".csv";
	_fileManager.create_file(_csvCheck);
	
	std::string _csvSubSuc = currTime + "_2hdmAlignedEdmAmmSuccess" + ".csv";
	_fileManager.create_file(_csvSubSuc);
	
	std::vector<std::string> _header;
	_header.push_back("tan(beta)");
	_header.push_back("M12");
	for(int lamb = 1; lamb <=7; ++lamb){
		std::string curren = std::to_string(lamb);
		std::string currLambda = "Lambda" + curren;
		_header.push_back(currLambda);
	}
	
	_header.push_back("aU");
	_header.push_back("aD");
	_header.push_back("aL");
	
	for(int h = 0; h < 3; ++h){
		std::string currMh = "mh_" + std::to_string(h + 1);
		_header.push_back(currMh);
	}
	
	_header.push_back("mh_C");
	
	_header.push_back("eEDM");
	_header.push_back("aAMM");

	_header.push_back("1-loop EDM contribution");
	_header.push_back("1-loop AMM contribution");

	_fileManager.add_line(_csvCheck, _header);
	_fileManager.add_line(_csvSubSuc, _header);
	
	_header.clear();
	
	int totalSuccEDM = 0;
	
	//int totalSuccRGE = 0;
	
	int totalSubSuccDeAmm = 0;
	
	for(int tems = 1; tems <= datapoints; ++tems){
		
		SM _sm;
	
		THDM _thdm(_sm);
	
		Base_generic gen;
		
		std::vector<complex<double>> outputLine;
		
		gen.beta = std::log(beta(rng));
		double monetwo = std::sqrt(m12(rng));
		gen.M12 = std::complex<double>(monetwo, 0.);
		gen.Lambda1 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda1 * lambda1);
		gen.Lambda2 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda2 * lambda2);
		gen.Lambda3 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda3 * lambda3);
		gen.Lambda4 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda4 * lambda4);
		std::complex<double> _randPhase_5 = std::exp(_i * lam5phase(rng));
		gen.Lambda5 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda5 * lambda5) * _randPhase_5;
		std::complex<double> _randPhase_6 = std::exp(_i * lam6phase(rng));
		gen.Lambda6 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda6 * lambda6) * _randPhase_6;
		std::complex<double> _randPhase_7 = std::exp(_i * lam7phase(rng));
		gen.Lambda7 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda7 * lambda7) * _randPhase_7;
		
		_thdm.set_param_gen(gen);
		
		
		std::complex<double> aU = upInter(rng) * std::exp(_i * upPhase(rng));
		std::complex<double> aD = downInter(rng) * std::exp(_i * downPhase(rng));
		std::complex<double> aL = leptonInter(rng) * std::exp(_i * leptonPhase(rng));
		
		_thdm.set_yukawa_aligned(aU, aD, aL);
		
		outputLine.push_back(std::tan(gen.beta));
		outputLine.push_back(gen.M12);
		outputLine.push_back(gen.Lambda1);
		outputLine.push_back(gen.Lambda2);
		outputLine.push_back(gen.Lambda3);
		outputLine.push_back(gen.Lambda4);
		outputLine.push_back(gen.Lambda5);
		outputLine.push_back(gen.Lambda6);
		outputLine.push_back(gen.Lambda7);
		
		outputLine.push_back(aU);
		outputLine.push_back(aD);
		outputLine.push_back(aL);
		
		vector<double> _higgsMasses = _thdm.get_higgs_treeLvl_masses();
		
		for(int higgs = 0; higgs <= 3; ++higgs){
			outputLine.push_back(_higgsMasses[higgs]);
		}
		
		
		if( (120 < _higgsMasses[0]) && (_higgsMasses[0] < 130)){
			std::cout << "Higgs mass between 120 and 130 GeV\n";
			
			EDM _edm(_thdm);

			if(_edm.is_initialized()){
				
				std::cout << "EDM type object initialised\n";
				
				std::cout << tems << " out of " << datapoints << " random tests done\n";
				
				if( (_thdm.is_perturbative() ) && (_thdm.is_unitary() ) && (_thdm.is_stable() ) ){
				
					totalSuccEDM += 1;
				
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
				
					outputLine.push_back(_eEDM);
					outputLine.push_back(_aAMM);
				
					outputLine.push_back(_eEDM_1loop);
					outputLine.push_back(_aAMM_1loop);
					
					_fileManager.add_line(_csvCheck, outputLine);
				
					if( (_eEDM < de) && (amu < _aAMM)){
						_fileManager.add_line(_csvSubSuc, outputLine);
						std::cout << "d_e and a_mu fulfilled input conditions\n";
					
						totalSubSuccDeAmm += 1;
					}
					
				}
			}
		}
		
		_higgsMasses.clear();
		
		outputLine.clear();
		
	}
	
	std::cout << "\n\n\nFinalised\n\n" << "Successfully initialised " << totalSuccEDM << " EDMs, out of which " << totalSubSuccDeAmm << " were found in interesting regions for d_e and a_mu.\n\n\n";
}



void genSearcher(double mintanbeta, double maxtanbeta, double minmonetwo, 
		double maxonetwo, double lambda1, double lambda2, double lambda3, double lambda4, double lambda5, double lam5minphase, double lam5maxphase, 
		double lambda6, double lam6minphase, double lam6maxphase, double lambda7, double lam7minphase, double lam7maxphase, 
		double updiagmin, double updiagmax, double updiagphasemin, double updiagphasemax, double upoffdiagmin, double upoffdiagmax, double upoffdiagphasemin, double upoffdiagphasemax, 
		double downdiagmin, double downdiagmax, double downdiagphasemin, double downdiagphasemax, double downoffdiagmin, double downoffdiagmax, double downoffdiagphasemin, double downoffdiagphasemax, 
		double leptondiagmin, double leptondiagmax, double leptondiagphasemin, double leptondiagphasemax, double leptonoffdiagmin, double leptonoffdiagmax, double leptonoffdiagphasemin, double leptonoffdiagphasemax, 
		double amu, double de, int datapoints){
	
	time_t rawtime;
	time (&rawtime);
	std::string curTime = ctime(&rawtime);
	std::string currTime = curTime.substr(0, curTime.size() - 1);
	currTime.erase(0,4);
	
	std::mt19937_64 rng;
    // initialize random number generator with time-dependent seed
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    std::uniform_real_distribution<double> unif(0, 1);
	std::uniform_real_distribution<double> sign(-1,1);
	std::uniform_real_distribution<double> beta(std::exp(std::atan(mintanbeta)), std::exp(std::atan(maxtanbeta)));
	std::uniform_real_distribution<double> m12(minmonetwo * minmonetwo, maxonetwo * maxonetwo);
	std::uniform_real_distribution<double> lam5phase(lam5minphase, lam5maxphase);
	std::uniform_real_distribution<double> lam6phase(lam5minphase, lam5maxphase);
	std::uniform_real_distribution<double> lam7phase(lam5minphase, lam5maxphase);
	std::uniform_real_distribution<double> upDiagInt(updiagmin, updiagmax);
	std::uniform_real_distribution<double> upDiagPhase(updiagphasemin, updiagphasemax);
	std::uniform_real_distribution<double> upOffDiagInt(upoffdiagmin, upoffdiagmax);
	std::uniform_real_distribution<double> upOffDiagPhase(upoffdiagphasemin, upoffdiagphasemax);
	std::uniform_real_distribution<double> downDiagInt(downdiagmin, downdiagmax);
	std::uniform_real_distribution<double> downDiagPhase(downdiagphasemin, downdiagphasemax);
	std::uniform_real_distribution<double> downOffDiagInt(downoffdiagmin, downoffdiagmax);
	std::uniform_real_distribution<double> downOffDiagPhase(downoffdiagphasemin, downoffdiagphasemax);
	std::uniform_real_distribution<double> leptonDiagInt(leptondiagmin, leptondiagmax);
	std::uniform_real_distribution<double> leptonDiagPhase(leptondiagphasemin, leptondiagphasemax);
	std::uniform_real_distribution<double> leptonOffDiagInt(leptonoffdiagmin, leptonoffdiagmax);
	std::uniform_real_distribution<double> leptonOffDiagPhase(leptonoffdiagphasemin, leptonoffdiagphasemax);
    // ready to generate random numbers
	
	FileSystem _fileManager;
	std::string _outDir = "output/data/";
	_fileManager.set_directory(_outDir);
	
	std::string _csvCheck = currTime + "_2hdmGenSuccessInit" + ".csv";
	_fileManager.create_file(_csvCheck);
	
	std::string _csvSubSuc = currTime + "_2hdmGenEdmAmmSuccess" + ".csv";
	_fileManager.create_file(_csvSubSuc);
	
	std::vector<std::string> _header;
	_header.push_back("tan(beta)");
	_header.push_back("M12");
	for(int lamb = 1; lamb <=7; ++lamb){
		std::string curren = std::to_string(lamb);
		std::string currLambda = "Lambda" + curren;
		_header.push_back(currLambda);
	}
	
	_header.push_back("rho^U diag");
	
	for(int row = 0; row < 3; ++row){
		for(int colum = 0; colum < 3; ++ colum){
		std::string currRho;
		currRho = "lambdaU_(" + std::to_string(row + 1) + "," + std::to_string(colum + 1) + ")";
		_header.push_back(currRho);
		}
	}
	
	for(int row = 0; row < 3; ++row){
		for(int colum = 0; colum < 3; ++ colum){
		std::string currRho;
		currRho = "lambdaD_(" + std::to_string(row + 1) + "," + std::to_string(colum + 1) + ")";
		_header.push_back(currRho);
		}
	}
	
	for(int row = 0; row < 3; ++row){
		for(int colum = 0; colum < 3; ++ colum){
		std::string currRho;
		currRho = "lambdaL_(" + std::to_string(row + 1) + "," + std::to_string(colum + 1) + ")";
		_header.push_back(currRho);
		}
	}
	
	for(int h = 0; h < 3; ++h){
		std::string currMh = "mh_" + std::to_string(h + 1);
		_header.push_back(currMh);
	}
	
	_header.push_back("mh_C");
	
	_header.push_back("eEDM");
	_header.push_back("aAMM");

	_header.push_back("1-loop EDM contribution");
	_header.push_back("1-loop AMM contribution");

	_fileManager.add_line(_csvCheck, _header);
	_fileManager.add_line(_csvSubSuc, _header);
	
	_header.clear();
	
	int totalSuccEDM = 0;
	
	int totalSubSuccDeAmm = 0;
	
	for(int tems = 1; tems <= datapoints; ++tems){
		
		SM _sm;
	
		THDM _thdm(_sm);
	
		Base_generic gen;
		
		std::vector<complex<double>> outputLine;
		
		gen.beta = std::log(beta(rng));
		double monetwo = std::sqrt(m12(rng));
		gen.M12 = std::complex<double>(monetwo, 0.);
		gen.Lambda1 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda1 * lambda1);
		gen.Lambda2 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda2 * lambda2);
		gen.Lambda3 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda3 * lambda3);
		gen.Lambda4 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda4 * lambda4);
		std::complex<double> _randPhase_5 = std::exp(_i * lam5phase(rng));
		gen.Lambda5 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda5 * lambda5) * _randPhase_5;
		std::complex<double> _randPhase_6 = std::exp(_i * lam6phase(rng));
		gen.Lambda6 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda6 * lambda6) * _randPhase_6;
		std::complex<double> _randPhase_7 = std::exp(_i * lam7phase(rng));
		gen.Lambda7 = signum(sign(rng)) * std::sqrt(unif(rng) * lambda7 * lambda7) * _randPhase_7;
		
		outputLine.push_back(std::tan(gen.beta));
		outputLine.push_back(gen.M12);
		outputLine.push_back(gen.Lambda1);
		outputLine.push_back(gen.Lambda2);
		outputLine.push_back(gen.Lambda3);
		outputLine.push_back(gen.Lambda4);
		outputLine.push_back(gen.Lambda5);
		outputLine.push_back(gen.Lambda6);
		outputLine.push_back(gen.Lambda7);
		
		Eigen::Matrix3cd rho_Ulambdas;
		Eigen::Matrix3cd rho_Dlambdas;
		Eigen::Matrix3cd rho_Llambdas;
		
		rho_Ulambdas.setZero();
		rho_Dlambdas.setZero();
		rho_Llambdas.setZero();
		
		double lambdaUdiag = upDiagInt(rng);
		double Udiagphase = upDiagPhase(rng);
		double lambdaDdiag = downDiagInt(rng);
		double Ddiagphase = downDiagPhase(rng);
		double lambdaLdiag = leptonDiagInt(rng);
		double Ldiagphase = leptonDiagPhase(rng);
		
		for(int rowpos = 0; rowpos < 3; ++rowpos){
			for(int colpos = 0; colpos < 3; ++colpos){
				if(rowpos == colpos){
					rho_Ulambdas(rowpos, colpos) = lambdaUdiag * std::exp(_i * Udiagphase);
				} else{
					rho_Ulambdas(rowpos, colpos) = upOffDiagInt(rng) * std::exp(_i * upOffDiagPhase(rng));
				}
				outputLine.push_back(rho_Ulambdas(rowpos, colpos));
			}
		}
		
		for(int rowpos = 0; rowpos < 3; ++rowpos){
			for(int colpos = 0; colpos < 3; ++colpos){
				if(rowpos == colpos){
					rho_Dlambdas(rowpos, colpos) = lambdaDdiag * std::exp(_i * Ddiagphase);
				} else{
					rho_Dlambdas(rowpos, colpos) = downOffDiagInt(rng) * std::exp(_i * downOffDiagPhase(rng));
				}
				outputLine.push_back(rho_Dlambdas(rowpos, colpos));
			}
		}
		
		for(int rowpos = 0; rowpos < 3; ++rowpos){
			for(int colpos = 0; colpos < 3; ++colpos){
				if(rowpos == colpos){
					rho_Llambdas(rowpos, colpos) = lambdaLdiag * std::exp(_i * Ldiagphase);
				} else{
					rho_Llambdas(rowpos, colpos) = leptonOffDiagInt(rng) * std::exp(_i * leptonOffDiagPhase(rng));
				}
				outputLine.push_back(rho_Llambdas(rowpos, colpos));
			}
		}
		
		Eigen::Matrix3cd rhoU;
		Eigen::Matrix3cd rhoD;
		Eigen::Matrix3cd rhoL;
		rhoU.setZero();
		rhoD.setZero();
		rhoL.setZero();
		
		double _v2 = _thdm.get_v2();
		
		FermionSector up = UP;
		FermionSector down = DOWN;
		FermionSector lepton = LEPTON;
		
		for(int rowpos = 0; rowpos < 3; ++rowpos){
			for(int colpos = 0; colpos < 3; ++colpos){
				if(rowpos == colpos){
					rhoU(rowpos, colpos) = rho_Ulambdas(rowpos, colpos) * _thdm.get_mf(up, rowpos + 1) * std::sqrt(2. / _v2);
					rhoD(rowpos, colpos) = rho_Dlambdas(rowpos, colpos) * _thdm.get_mf(down, rowpos + 1) * std::sqrt(2. / _v2);
					rhoL(rowpos, colpos) = rho_Llambdas(rowpos, colpos) * _thdm.get_mf(lepton, rowpos + 1) * std::sqrt(2. / _v2);
				} else{
					rhoU(rowpos, colpos) = rho_Ulambdas(rowpos, colpos) * std::sqrt(_thdm.get_mf(up, rowpos + 1) * _thdm.get_mf(up, colpos + 1) * 2. / _v2);
					rhoD(rowpos, colpos) = rho_Dlambdas(rowpos, colpos) * std::sqrt(_thdm.get_mf(down, rowpos + 1) * _thdm.get_mf(down, colpos + 1) * 2. / _v2);
					rhoL(rowpos, colpos) = rho_Llambdas(rowpos, colpos) * std::sqrt(_thdm.get_mf(lepton, rowpos + 1) * _thdm.get_mf(lepton, colpos + 1) * 2. / _v2);
				}
			}
		}
		
		_thdm.set_param_gen(gen);
		
		vector<double> _higgsMasses = _thdm.get_higgs_treeLvl_masses();
		
		for(int higgs = 0; higgs <= 3; ++higgs){
			outputLine.push_back(_higgsMasses[higgs]);
		}
		
		
		if( (120 < _higgsMasses[0]) && (_higgsMasses[0] < 130)){
			std::cout << "Higgs mass between 120 and 130 GeV\n";
			
			EDM _edm(_thdm);

			if(_edm.is_initialized()){
				
				std::cout << "EDM type object initialised\n";
				
				std::cout << tems << " out of " << datapoints << " random tests done\n";
				
				if( (_thdm.is_perturbative() ) && (_thdm.is_unitary() ) && (_thdm.is_stable() ) ){
				
					totalSuccEDM += 1;
				
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
				
					outputLine.push_back(_eEDM);
					outputLine.push_back(_aAMM);
				
					outputLine.push_back(_eEDM_1loop);
					outputLine.push_back(_aAMM_1loop);
					
					_fileManager.add_line(_csvCheck, outputLine);
				
					if( (_eEDM < de) && (amu < _aAMM)){
						_fileManager.add_line(_csvSubSuc, outputLine);
						std::cout << "d_e and a_mu fulfilled input conditions\n";
					
						totalSubSuccDeAmm += 1;
					}
					
				}
			}
		}
		
		_higgsMasses.clear();
		
		outputLine.clear();
		
	}
	
	std::cout << "\n\n\nFinalised\n\n" << "Successfully initialised " << totalSuccEDM << " EDMs, out of which " << totalSubSuccDeAmm << " were found in interesting regions for d_e and a_mu.\n\n\n";
}



int main(){
	
	double betamin = 0.;
	double betamax = 0.;
	double m12min = 0.;
	double m12max = 0.;
	double lambda1max = 0.;
	double lambda2max = 0.;
	double lambda3max = 0.;
	double lambda4max = 0.;
	double lambda5max = 0.;
	double lambda5phasemin = 0.;
	double lambda5phasemax = 0.;
	double lambda6max = 0.;
	double lambda6phasemin = 0.;
	double lambda6phasemax = 0.;
	double lambda7max = 0.;
	double lambda7phasemin = 0.;
	double lambda7phasemax = 0.;
	double edmmax = 0.;
	double ammmin = 0.;
	int datanr = 0;
	
	
	std::cout << "\n\nWelcome to the 2hdmSearcher program. This program allows you to input parameter regions and search them randomly for 2hdm candidate points based on their a_mu and d_e contributions. Note that there is little protection for incorrect inputs, so beware.\n\n";
	std::cout << "Would you like to test a hard or soft Z_2-symmetric 2HDM? (Y/N)\n";
	std::string typeresp;
	std::getline (std::cin,typeresp);
	
	if((typeresp == "Y") or (typeresp == "y") or (typeresp == "yes") or (typeresp == "Yes")){
		std::cout << "\n" << "Lower bound for tan(beta): ";
		std::cin >> betamin;
		std::cout << "Upper bound for tan(beta): ";
		std::cin >> betamax;
		std::cout << "Lower bound for |m_12^2|: ";
		std::cin >> m12min;
		std::cout << "Upper bound for |m_12^2|: ";
		std::cin >> m12max;
		std::cout << "Upper bound for |lambda_1|: ";
		std::cin >> lambda1max;
		std::cout << "Upper bound for |lambda_2|: ";
		std::cin >> lambda2max;
		std::cout << "Upper bound for |lambda_3|: ";
		std::cin >> lambda3max;
		std::cout << "Upper bound for |lambda_4|: ";
		std::cin >> lambda4max;
		std::cout << "Upper bound for |lambda_5|: ";
		std::cin >> lambda5max;
		std::cout << "Lower bound for phase of lambda_5: ";
		std::cin >> lambda5phasemin;
		std::cout << "Upper bound for phase of lambda_5: ";
		std::cin >> lambda5phasemax;
		//std::cout << "\n" << "Upper bound for |lambda_6|: ";
		//std::cin >> lambda6max;
		//std::cout << "\n" << "Lower bound for phase of lambda_6: ";
		//std::cin >> lambda6phasemin;
		//std::cout << "\n" << "Upper bound for phase of lambda_6: ";
		//std::cin >> lambda6phasemax;
		//std::cout << "\n" << "Upper bound for |lambda_7|: ";
		//std::cin >> lambda7max;
		//std::cout << "\n" << "Lower bound for phase of lambda_7: ";
		//std::cin >> lambda7phasemin;
		//std::cout << "\n" << "Upper bound for phase of lambda_7: ";
		//std::cin >> lambda7phasemax;
		std::cout << "\n" << "Upper bound for d_e: ";
		std::cin >> edmmax;
		std::cout << "\n" << "Lower bound for a_mu: ";
		std::cin >> ammmin;
		
		std::cout << "\n\n" << "How many random datapoints would you like to check? \n";
		std::cin >> datanr;

		std::string typetest;
		std::cout << "\n\nNow you may decide what type of 2HDM you want to explore. \n \nDo you wish to test Type I, II, Y, or X? (input as \"I\", \"II\", \"Y\", or \"X\", without the quotation marks) \n";
		
		std::getline (std::cin,typetest);
		std::getline (std::cin,typetest);
		
		if((typetest == "I") or (typetest == "1")){
			typeSearcher("I", betamin, betamax, m12min, m12max, lambda1max, lambda2max, lambda3max, lambda4max, 
	lambda5max, lambda5phasemin, lambda5phasemax, ammmin, edmmax, datanr);
		} else if((typetest == "Y") or (typetest == "III")){
			typeSearcher("Y", betamin, betamax, m12min, m12max, lambda1max, lambda2max, lambda3max, lambda4max, 
	lambda5max, lambda5phasemin, lambda5phasemax, ammmin, edmmax, datanr);
		} else if((typetest == "II") or (typetest == "2")){
			typeSearcher("II", betamin, betamax, m12min, m12max, lambda1max, lambda2max, lambda3max, lambda4max, 
	lambda5max, lambda5phasemin, lambda5phasemax, ammmin, edmmax, datanr);
		} else if((typetest == "X") or (typetest == "IV")){
			typeSearcher("X", betamin, betamax, m12min, m12max, lambda1max, lambda2max, lambda3max, lambda4max, 
	lambda5max, lambda5phasemin, lambda5phasemax, ammmin, edmmax, datanr);
		}
	} else{
		std::cout << "\n\n" << "Would you like to test an alligned 2HDM? (Y/N)\n";
		std::string allignresp;
		std::getline (std::cin,allignresp);
	
		if((allignresp == "Y") or (allignresp == "y") or (allignresp == "yes") or (allignresp == "Yes")){
			std::cout << "\n" << "Lower bound for tan(beta): ";
			std::cin >> betamin;
			std::cout << "Upper bound for tan(beta): ";
			std::cin >> betamax;
			std::cout << "Lower bound for |m_12^2|: ";
			std::cin >> m12min;
			std::cout << "Upper bound for |m_12^2|: ";
			std::cin >> m12max;
			std::cout << "Upper bound for |lambda_1|: ";
			std::cin >> lambda1max;
			std::cout << "Upper bound for |lambda_2|: ";
			std::cin >> lambda2max;
			std::cout << "Upper bound for |lambda_3|: ";
			std::cin >> lambda3max;
			std::cout << "Upper bound for |lambda_4|: ";
			std::cin >> lambda4max;
			std::cout << "Upper bound for |lambda_5|: ";
			std::cin >> lambda5max;
			std::cout << "Lower bound for phase of lambda_5: ";
			std::cin >> lambda5phasemin;
			std::cout << "Upper bound for phase of lambda_5: ";
			std::cin >> lambda5phasemax;
			std::cout << "Upper bound for |lambda_6|: ";
			std::cin >> lambda6max;
			std::cout << "Lower bound for phase of lambda_6: ";
			std::cin >> lambda6phasemin;
			std::cout << "Upper bound for phase of lambda_6: ";
			std::cin >> lambda6phasemax;
			std::cout << "Upper bound for |lambda_7|: ";
			std::cin >> lambda7max;
			std::cout << "Lower bound for phase of lambda_7: ";
			std::cin >> lambda7phasemin;
			std::cout << "Upper bound for phase of lambda_7: ";
			std::cin >> lambda7phasemax;
			std::cout << "Upper bound for d_e: ";
			std::cin >> edmmax;
			std::cout << "Lower bound for a_mu: ";
			std::cin >>ammmin;
			
			double leptalmin = 0.;
			double leptalmax = 0.;
			double upalmin = 0.;
			double upalmax = 0.;
			double downalmin = 0.;
			double downalmax = 0.;
			double leptphasemin = 0.;
			double leptphasemax = 0.;
			double upphasemin = 0.;
			double upphasemax = 0.;
			double downphasemin = 0.;
			double downphasemax = 0.;
			
		std::cout << "\n\nNow you will choose what alignment region you wish to explore.\n\n";
			
			std::cout << "Lower bound for a^U: ";
			std::cin >> upalmin;
			std::cout << "Upper bound for a^U: ";
			std::cin >> upalmax;
			std::cout << "Lower bound for phase of a^U: ";
			std::cin >> upphasemin;
			std::cout << "Upper bound for phase of a^U: ";
			std::cin >> upphasemax;
			std::cout << "Lower bound for a^D: ";
			std::cin >> downalmin;
			std::cout << "Upper bound for a^D: ";
			std::cin >> downalmax;
			std::cout << "Lower bound for phase of a^D: ";
			std::cin >> downphasemin;
			std::cout << "Upper bound for phase of a^D: ";
			std::cin >> downphasemax;
			std::cout << "Lower bound for a^L: ";
			std::cin >> leptalmin;
			std::cout << "Upper bound for a^L: ";
			std::cin >> leptalmax;
			std::cout << "Lower bound for phase of a^L: ";
			std::cin >> leptphasemin;
			std::cout << "Upper bound for phase of a^L: ";
			std::cin >> leptphasemax;
			
			std::cout << "\nHow many datapoints would you like to check?\n";
			std::cin >> datanr;
			
			std::cout <<"\n\n";
			
			alignSearcher(betamin, betamax, m12min, 
		m12max, lambda1max, lambda2max, lambda3max, lambda4max, lambda5max, lambda5phasemin, lambda5phasemax, 
		lambda6max, lambda6phasemin, lambda6phasemax, lambda7max, lambda7phasemin, lambda7phasemax, 
		upalmin, upalmax, upphasemin, upphasemax, downalmin, downalmax, downphasemin, downphasemax, 
		leptalmin, leptalmax, leptphasemin, leptphasemax, ammmin, edmmax, datanr);
		} else{
			
			std::cout << "You have chosen to test a fully general 2HDM parameter space.\n\n";
			
			std::cout << "\n" << "Lower bound for tan(beta): ";
			std::cin >> betamin;
			std::cout << "Upper bound for tan(beta): ";
			std::cin >> betamax;
			std::cout << "Lower bound for |m_12^2|: ";
			std::cin >> m12min;
			std::cout << "Upper bound for |m_12^2|: ";
			std::cin >> m12max;
			std::cout << "Upper bound for |lambda_1|: ";
			std::cin >> lambda1max;
			std::cout << "Upper bound for |lambda_2|: ";
			std::cin >> lambda2max;
			std::cout << "Upper bound for |lambda_3|: ";
			std::cin >> lambda3max;
			std::cout << "Upper bound for |lambda_4|: ";
			std::cin >> lambda4max;
			std::cout << "Upper bound for |lambda_5|: ";
			std::cin >> lambda5max;
			std::cout << "Lower bound for phase of lambda_5: ";
			std::cin >> lambda5phasemin;
			std::cout << "Upper bound for phase of lambda_5: ";
			std::cin >> lambda5phasemax;
			std::cout << "Upper bound for |lambda_6|: ";
			std::cin >> lambda6max;
			std::cout << "Lower bound for phase of lambda_6: ";
			std::cin >> lambda6phasemin;
			std::cout << "Upper bound for phase of lambda_6: ";
			std::cin >> lambda6phasemax;
			std::cout << "Upper bound for |lambda_7|: ";
			std::cin >> lambda7max;
			std::cout << "Lower bound for phase of lambda_7: ";
			std::cin >> lambda7phasemin;
			std::cout << "Upper bound for phase of lambda_7: ";
			std::cin >> lambda7phasemax;
			std::cout << "Upper bound for d_e: ";
			std::cin >> edmmax;
			std::cout << "Lower bound for a_mu: ";
			std::cin >>ammmin;
			
			double updiagmin = 0.;
			double updiagmax = 0.;
			double updiagphasemin = 0.;
			double updiagphasemax = 0.;
			double upoffdiagmin = 0.;
			double upoffdiagmax = 0.;
			double upoffdiagphasemin = 0.;
			double upoffdiagphasemax = 0.;
			double downdiagmin = 0.;
			double downdiagmax = 0.;
			double downdiagphasemin = 0.;
			double downdiagphasemax = 0.;
			double downoffdiagmin = 0.;
			double downoffdiagmax = 0.;
			double downoffdiagphasemin = 0.;
			double downoffdiagphasemax = 0.;
			double leptondiagmin = 0.;
			double leptondiagmax = 0.;
			double leptondiagphasemin = 0.;
			double leptondiagphasemax = 0.;
			double leptonoffdiagmin = 0.;
			double leptonoffdiagmax = 0.;
			double leptonoffdiagphasemin = 0.;
			double leptonoffdiagphasemax = 0.;
			
			std::cout << "\n\nNext, you will set the rho^F parameter regions, using the Cheng-Sher ansatz. Note that all diagonal lambdas will be identical for a given parameter point.\n\n";
			
			std::cout << "Lower bound for up-sector diagonal: ";
			std::cin >> updiagmin;
			std::cout << "Upper bound for up-sector diagonal: ";
			std::cin >> updiagmax;
			std::cout << "Lower bound for up-sector diagonal phase: ";
			std::cin >> updiagphasemin;
			std::cout << "Upper bound for up-sector diagonal phase: ";
			std::cin >> updiagphasemax;
			std::cout << "Lower bound for up-sector off-diagonal: ";
			std::cin >> upoffdiagmin;
			std::cout << "Upper bound for up-sector off-diagonal: ";
			std::cin >> upoffdiagmax;
			std::cout << "Lower bound for up-sector off-diagonal phase: ";
			std::cin >> upoffdiagphasemin;
			std::cout << "Upper bound for up-sector off-diagonal phase: ";
			std::cin >> upoffdiagphasemax;
			std::cout << "Lower bound for down-sector diagonal: ";
			std::cin >> downdiagmin;
			std::cout << "Upper bound for down-sector diagonal: ";
			std::cin >> downdiagmax;
			std::cout << "Lower bound for down-sector diagonal phase: ";
			std::cin >> downdiagphasemin;
			std::cout << "Upper bound for down-sector diagonal phase: ";
			std::cin >> downdiagphasemax;
			std::cout << "Lower bound for down-sector off-diagonal: ";
			std::cin >> downoffdiagmin;
			std::cout << "Upper bound for down-sector off-diagonal: ";
			std::cin >> downoffdiagmax;
			std::cout << "Lower bound for down-sector off-diagonal phase: ";
			std::cin >> downoffdiagphasemin;
			std::cout << "Upper bound for down-sector off-diagonal phase: ";
			std::cin >> downoffdiagphasemax;
			std::cout << "Lower bound for lepton-sector diagonal: ";
			std::cin >> leptondiagmin;
			std::cout << "Upper bound for lepton-sector diagonal: ";
			std::cin >> leptondiagmax;
			std::cout << "Lower bound for lepton-sector diagonal phase: ";
			std::cin >> leptondiagphasemin;
			std::cout << "Upper bound for lepton-sector diagonal phase: ";
			std::cin >> leptondiagphasemax;
			std::cout << "Lower bound for lepton-sector off-diagonal: ";
			std::cin >> leptonoffdiagmin;
			std::cout << "Upper bound for lepton-sector off-diagonal: ";
			std::cin >> leptonoffdiagmax;
			std::cout << "Lower bound for lepton-sector off-diagonal phase: ";
			std::cin >> leptonoffdiagphasemin;
			std::cout << "Upper bound for lepton-sector off-diagonal phase: ";
			std::cin >> leptonoffdiagphasemax;
			
			std::cout << "\n\nHow many datapoints would you like to test?\n";
			std::cin >> datanr;
			std::cout << "\n\n";
			
			genSearcher(betamin, betamax, m12min, 
		m12max, lambda1max, lambda2max, lambda3max, lambda4max, lambda5max, lambda5phasemin, lambda5phasemax, 
		lambda6max, lambda6phasemin, lambda6phasemax, lambda7max, lambda7phasemin, lambda7phasemax, 
		updiagmin, updiagmax, updiagphasemin, updiagphasemax, upoffdiagmin, upoffdiagmax, upoffdiagphasemin, upoffdiagphasemax, 
		downdiagmin, downdiagmax, downdiagphasemin, downdiagphasemax, downoffdiagmin, downoffdiagmax, downoffdiagphasemin, downoffdiagphasemax, 
		leptondiagmin, leptondiagmax, leptondiagphasemin, leptondiagphasemax, leptonoffdiagmin, leptonoffdiagmax, leptonoffdiagphasemin, leptonoffdiagphasemax, 
		ammmin, edmmax, datanr);
			
		}
	}
	
	
	std::cout << "\n\nProgram finished.\n\n";
	
}




