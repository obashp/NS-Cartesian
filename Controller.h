#include "ColDriver.h"

class Controller
{
	Geometry *G;
	MomDriver *M;
	PresDriver *P;
	ScalDriver *E;

	double rho, mu, Cp, gx, gy, Beta, Pr, dt;
	double *Lambda;
	double *Tolerance;
	double *Urf;
	double Coupling;
	long int MaxIterations;
	long int *Inners;
	unsigned short Time_Scheme;
	long int Ntime, DataIntervals, TimeStep;
	unsigned short *eq;
	unsigned short Neq;
	unsigned short RestartFlag, buoyancy;
	unsigned short StoreResiduals;
	unsigned short SolverCode;
	fstream f;

	static long int iOuter;

	bool PrintResiduals();
	void Read_InputFile();
public:
	Controller();
	~Controller();
	inline static void IncOuter(){iOuter = iOuter+1; Driver::IncOuter();}
	inline static long int GetOuter(){return iOuter;}
	void Initialize();
	void Run();
};


