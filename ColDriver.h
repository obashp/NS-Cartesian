#include "Grid.h"
#include "Solver.h"
#pragma once
using namespace std;

#ifndef COLDRIVER_H
#define COLDRIVER_H


class Driver;
class MomDriver;
class PresDriver;

enum dir
{X = 0, Y = 1};

/*
	Driver class is an abstract class, which contains variables required for
	driving a generic convection-diffusion equation.
	The Navier-Stokes equations are formulated as derived classes of the Driver
*/
class Driver
{
protected:
	Geometry *G;
	unsigned short nDiag;	//No. of diagonals in A
	unsigned short nvars;
	double **A, **S;	//S takes into consideration the number of variables
	double **phi, **phi_, **phi__;	//Solution vector, and previous time solutions in case of unsteady problems
	int *NbrIdx;	//Neighbour indices of P- can be a 1D array for Structured Grid
	double Cn, Ce, Cp, D;	//Convective and Diffusive Flux Coefficients
	double *Lambda, *Urf;	//Blending factor for Deferred Correction, Underrelaxation factor for diagonal dominance
	double Tol;
	long int inners;
	//timestep, physical time, time instant  
	static unsigned short restart;
	static double dt, time;
	static unsigned long int ntime, itime;
	static unsigned short Time_Scheme;	// 0 - Euler Implicit, 1- BDF2, 2- Steady

	Solver *Sol;
	virtual void SetBC() = 0;
	virtual void ComputeConDiff() = 0;
	virtual void SetEquation(int varc) = 0;
	void StoreMat(string fn);
	void StoreS(string fn, int varc);
public:
	Driver();
	~Driver();
	void SetGeometry(Geometry *in_G) {G = in_G;}
	virtual void Initialize(Geometry *inG, int in_nvars, int in_nDiag = 5);
	inline void SetLambda(double iLambda, int varc){Lambda[varc] = iLambda;}
	inline void SetUrf(double iLambda, int varc){Urf[varc] = iLambda;}
	inline void SetTolerance(double Tolerance){Tol = Tolerance;}
	inline void SetInners(long int Inners){inners = Inners;}
	inline static void Setrestart(unsigned short res){restart = res;}
	inline static unsigned short Getrestart(){return restart;}
	inline static void SetTStep(double idt){dt = idt;}
	inline static void SetNtime(unsigned long int intime){ntime = intime;}
	inline static void IncTime(){time = time+dt; itime++;}
	inline static unsigned long int GetiTime(){return itime;}
	inline static double GetTimeStep(){return dt;}
	inline static double GetTime(){return time;}
	inline static void SetTimeScheme(int iTScheme){Time_Scheme = iTScheme;}
	inline static unsigned short GetTimeScheme(){return Time_Scheme;}

	const double* GetPhi(int in_var){return phi[in_var];}
	virtual void Solve() = 0;
	virtual void Init() = 0;
	virtual void Store() = 0;
	
};

class MomDriver:public Driver
{
	unsigned short bodyforcing;
	// Body forcing variable- obtained from corresponding class
	const double *vF;
	// pressure vector-obtained from pressure driver class
	// mdotE, mdotN - mass flux vectors- obtained from pressure driver class
	const double *p, *mdotE, *mdotN;	
	//Diagonal entries for applying SIMPLE Algorithm
	double **Ap, **pgrad;
	//density, viscosity
	double rho, mu;
	//gravity, compressibility factor
	double gx, gy, beta;
	double MomRes[2];

	enum Vel
	{U = 0, V = 1};

	void SetBC();
	void SetMomBC(int varc);
	void ComputePGrad();
	void ComputeSource();
	void ComputeConDiff();
	void ComputeTime();
	void SetEquation(int varc);
public:
	MomDriver(double irho, double imu);
	MomDriver(double irho, double imu, double igx, double igy, double ibeta); 
	void Initialize(Geometry *inG, int nDim = 2, int nDiag = 5);
	inline void SetPressure(const double *in_p){p = in_p;}
	inline void SetmdotE(const double *in_mE){mdotE = in_mE;}
	inline void SetmdotN(const double *in_mN){mdotN = in_mN;}
	inline void SetFVar(const double *in_vF){vF = in_vF;}
	inline void SetTimeStep(double idt){dt = idt;}
	inline double* GetDiag(int varc){return Ap[varc];}
	inline double* GetPGrad(int varc){return pgrad[varc];}
	inline double* GetVel(int varc){return phi[varc];}
	inline double GetRes(int varc){return MomRes[varc];}
	inline void Init(){SetBC();}
	void Solve();
	void Store();
	
};


//A class for ensuring mass conservation, mass flux related boundary conditions are applied here
class PresDriver:public Driver
{
	//Additional Computables from pressure correction
 	double *p, *mdotE, *mdotN, *uf, *vf;
 	//Pressure gradient 
 	const double *pgradx, *pgrady, *Apu, *Apv, *u, *v;
 	// Density, reference pressure
 	double rho, pcref;
 	// Mass Correction and Total mass flow rates for Checking conservation
 	double MdotC, Mdot;
 	double Pres;

 	void SetBC();
 	void SetPresBC();
 	void SetEquation(int varc = 0);
 	void ComputeConDiff();
 	void ComputeSource();
 	void UpdateBoundaryPres(int varc);
 	void Corrections();
 public:
 	PresDriver(double irho);
 	~PresDriver();
 	void Initialize(Geometry *inG, int nvar = 1, int nDiag = 5);
 	inline void SetPGrad(const double* dpdl, int varc){if(varc == X) pgradx = dpdl; if(varc == Y) pgrady = dpdl;}
 	inline void SetAp(const double *Ap, int varc){if(varc == X) Apu = Ap; if(varc == Y) Apv = Ap;}
 	inline void SetVel(const double *iV, int varc){if(varc == X) u = iV; if(varc == Y) v = iV;}
 	inline void SetVelC(double *iV, int varc){if(varc == X) uf = iV; if(varc == Y) vf = iV;}
 	inline double* GetmdotE(){return mdotE;}
 	inline double* GetmdotN(){return mdotN;}
 	inline double* GetPressure(){return p;}
 	inline double GetRes(){return Pres;}
 	void Init();
 	void Solve();
 	void Store();
};

class ScalDriver:public Driver
{
	const double *mdotN, *mdotE;
	//density, viscosity, Prandtl Number
	double rho, mu, Pr, Cp;
	//Scalar Residual
	double Temp;

	enum Scal
	{T = 0};
	void SetBC();
	void SetScalBC();
	void ComputeSource();
	void ComputeConDiff();
	void ComputeTime();
	void SetEquation(int varc = T);
public:
	ScalDriver(double rho, double mu, double Pr, double Cp = 1.0);
	~ScalDriver();
	void Initialize(Geometry *inG, int nvar = 1, int nDiag = 5);
	inline void SetmdotE(const double *in_mE){mdotE = in_mE;}
	inline void SetmdotN(const double *in_mN){mdotN = in_mN;}
	void Init();
	inline double GetRes(){return Temp;}
	void Solve();
 	void Store();
};





#endif