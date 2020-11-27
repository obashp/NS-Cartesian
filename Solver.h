#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <string.h>
#include <sstream>

#define TOL 1e-12

using namespace std;

#ifndef SOLVER_H
#define SOLVER_H

void Allocate1Dvect(double **V, long int N, double def = 0.0);
void DeAllocate1Dvect(double **V, long int N);
void Reset1Dvect(double **V, long int N);
void Allocate2Dmat(double ***F, long int Nx, long int Ny, double def = 0.0);
void DeAllocate2Dmat(double ***F, long int Nx, long int Ny);

struct CSRmat
{
	long int Ncolumns;
	long int* columnoffset;
	long int* rows;
	double* values;

	CSRmat();
	CSRmat(long int iNrows);
	~CSRmat();
	void input(long int *icolumns, long int *irows, double *ivals);
	void Storemat(string fname);
};


class Solver
{
	double orf;	//Over-relaxation factor for Gauss-Seidel
	unsigned short Scode;	//Solver code
	unsigned short Pcode;	//Preconditioner code
	static unsigned int m;
	double Tolerance;

	const double *const *A;
	const double *b;
	const long int *Map21;	//Column dependence on Rows (Later to be replaced by CSRmat)
	const int *Nbr;	//Neighbour indices of index P, used for values of P
	unsigned short nNbr;	//Size of Nbr
	long int lNx, uNx;	//Indices to navigate through Map21
	long int N;	//Total Rows in Matrix

	double **LU, alpha;	//Parameters for SIP Preconditioner
	double *en, *en_, *Rn;	//Residual and error
	double *phi, *phi_;
	double ResNorm;
	double FinRes;
	int inneriters;

	// Variables for BiCGStab Solver
	double *v,*w, *p_, *s_, *r_;
	double rho1,rho2, Alpha, Beta, Omega;

	double *r_p, *p, *Ap, *r1, *s, *As, *s2;
	double a1, a2, a3;

	bool IsConverged(double *V, int N);
	void ComputeResidual();
	void SolveLU(double *x, double *b);
	void ILU();
	int GSSORSolver();
	int BiCGStabSolver();
	int PBiCGStabSolver();
public:
	Solver(){}
	~Solver();
	void Initialize(long int *Map, long int l,long int u, long int N);
	inline void SetnInner(int inner){inneriters = inner;}
	inline void SetA(const double *const *iA, const int *iNbr, unsigned short inNbr){A = iA; Nbr = iNbr; nNbr = inNbr;}
	inline void Setb(const double *ib){b = ib;}
	inline void Setphi(double *iphi){phi = iphi;}
	inline void SetTol(double iTol){Tolerance = iTol;}
	inline static void im(){m = m+1;}
	inline static unsigned int Getm(){return m;}
	inline void SetScode(unsigned short S){Scode = S;}
	inline void SetPcode(unsigned short P){Pcode = P;}
	inline double GetRes(){return ResNorm;}
	static double L1Norm(double *V, long int N);
	static double L2Norm(double *V, long int N);
	static double LinfNorm(double *V, long int N);
	static double dotp(double *V1, double *V2, long int N);
	static void copyvec(double *S, double *T, long int N);
	static void vecsum(double *C, double a,const double *A, double b, const double *B, long int N);
	static void vecscale(double *B, double a, double*A, long int N);
	static double sumvec(double *A, long int N);
	static void matvecmul(double *y, const double *const *A, const double *x, long int N, int nNbr, const int *Nbr);
	void Solve();

};

#endif
