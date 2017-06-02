#pragma once
#include "basefunctions.h"
#include "hdf5wrapper.h"
#include <gsl/gsl_integration.h>
typedef	std::vector<double> vec;



struct DynamicModelRK33{
	virtual ~DynamicModelRK33();
	virtual int f(double t, double *y, double *rhs) = 0;
	virtual void setManipulated(double *Tjkt) = 0; 
	// virtual void adaptSys() = 0; 

	double A2 = -5.0/9.0; 
	double A3= -153.0/128.0; 
	double B1 = 1./3.; 
	double B2 = 15./16.; 
	double B3 = 8./15.0;
	int step(double *k1, double *y, 
				double *rhs, double t, double h, uint ny, uint i0 = 0);

	
	vec y0;
	uint ny;
	vec tvector;
	double *pk1 = NULL;
	double *pyiter = NULL;
	double *prhs = NULL;
	double *pdsim = NULL;
	double *pmvspan = NULL;
	int setSimVariables(uint ny, vec y0, vec tsp, uint nu, double *mvsp);
	int sim(double *k1, double *y, 
				double *rhs, uint ny, vec const &tspan, 
				double *mvspan, double *dsim);


};



struct indx{
	uint C = 0;
	uint T = 1;
	std::vector<uint> mu = {2, 3, 4, 5};
};




struct CrystRawllignsModel : public DynamicModelRK33{

	CrystRawllignsModel();
	~CrystRawllignsModel();

	indx ii;
	double kb = 285.0; // 1/(s mu m^3)
	double b = 1.45;
	double EbbyR = 7517.0; // K
	double kg = 1.44e8; //mu m/s
	double g = 1.5;
	double EgbyR = 4859.0; //K
	double A = 0.25; //m^2
	double cp = 3.8; //kJ/K / kg
	double rho = 2.66e-12; // g/mu m^3
	double U = 1800.0 /3600.0; //kJ/m2 h K -> kJ/m2 s K
	double DeltaH = 44.5; //kJ/kg
	double M = 27.0; //kg
	double kv = 1.5;

	double nucleationRate(double T, double S, double mu3);
	double growthRate(double T, double S);
	double saturationComp(double T);
	double metastableComp(double T);
    static double PSDInititalFunction(double l, void * params);
    gsl_function F;    
    uint NpivActive;
	int NiICCalculation(std::vector<double> const &liBase, 
			uint Npivot_, std::vector<double> &NiIC, std::vector<double> &piIC);
 	vec referenceInitialState();
	void setManipulated(double *Tjkt); 

	 // some states to keep track:
	 double T, S, mu0, mu3, Cs, G, B0;

	 // Manipulated variable:
	 double Tjkt;
	
	 // Overriden RHS Function
	int f(double t, double *y, double *rhs);

	uint Nmu = 4;
	basefunc::MeasureTime tmr;

};

struct ParamEst
{
	ParamEst(){}
	~ParamEst(){}

	vec x0;
	vec xOtm, gOtm;
	double FOtm;
	vec dataEXP;
};