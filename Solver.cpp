#include "Solver.h"

using namespace std;

unsigned int Solver::m = 0;

void Allocate2Dmat(double ***F,long int Nx,long int Ny, double def)
{
	*F = new double*[Nx];
	for(long int I = 0; I < Nx; I++)
	{
		(*F)[I] = new double[Ny];
		for(long int J = 0; J < Ny; J++)
			(*F)[I][J] = def;
	}

}

void DeAllocate2Dmat(double ***F,long int Nx,long int Ny)
{
	for(int I = 0; I < Nx; I++)
	{
		delete[] (*F)[I];
		(*F)[I] = NULL;
	}

	delete[] *F;
	*F = NULL;
}

void Allocate1Dvect(double **V,long int N, double def)
{
	*V = new double[N];
	for(int I = 0; I < N; I++)
		(*V)[I] = def;
}

void DeAllocate1Dvect(double **V,long int N)
{
	delete[] *V;
	*V = NULL;
}

void Reset1Dvect(double **V,long int N)
{
	for(int I = 0 ; I < N; I++)
		(*V)[I] = 0.0;
}

CSRmat::CSRmat()
{
	Ncolumns = 0;
	columnoffset = NULL;
	rows = NULL; values = NULL;
}

CSRmat::CSRmat(long int iNcolumns) : Ncolumns(iNcolumns)
{
	columnoffset = new long int[Ncolumns+1];
	rows = NULL;
	values = NULL;
}

CSRmat::~CSRmat()
{
	if(values)	delete[] values;
	if(rows)	delete[] rows;
	if(columnoffset)	delete[] columnoffset;

	values = NULL; rows = NULL, columnoffset = NULL;
	Ncolumns = 0;
}

void CSRmat::input(long int *icolumns, long int *irows, double *ivals)
{
	for(int i = 0; i <= Ncolumns; i++)
		columnoffset[i] = icolumns[i];

	rows = new long int[columnoffset[Ncolumns]];
	values = new double[columnoffset[Ncolumns]];

	for(int i = 0; i < Ncolumns; i++)
		for(int j = columnoffset[i]; j < columnoffset[i+1]; j++)
		{
			rows[j] = irows[j];
			values[j] = ivals[j];
		}
}

void CSRmat::Storemat(string fname)
{
	fstream f;
	f.open(fname, ios::out);
	for(int i = 0 ; i < Ncolumns; i++)
	{
		for(int j = columnoffset[i]; j < columnoffset[i+1]; j++)
			f<<columnoffset[i]<<"	"<<rows[j]<<"	"<<values[j]<<endl;
		f<<endl;
	}
	f.close();
}

double Solver::L1Norm(double *V,long int N)
{
	double S = 0.0;
	for(int i = 0; i < N; i++)
		S += fabs(V[i]);
	return S; 
}

double Solver::L2Norm(double *V,long int N)
{
	return sqrt(dotp(V,V,N));
}

double Solver::LinfNorm(double *V, long int N)
{
	double S = fabs(V[0]);
	for(int i = 1; i < N; i++)
		if(fabs(V[i]) > S)
			S = fabs(V[i]);
	return S;
}

double Solver::sumvec(double *V,long int N)
{
	double Sum = 0.0;
	for(int I = 0; I < N; I++)
		Sum += (V[I]);
	return Sum;
}

double Solver::dotp(double *V1, double *V2,long int N)
{
	double S = 0.0;
	for(int i = 0; i < N; i++)
		S += V1[i]*V2[i];
	return S;
}

void Solver::copyvec(double *S, double *T,long int N)
{
	if(T == NULL) Allocate1Dvect(&T, N);

	for(int i = 0; i < N; i++)
		T[i] = S[i];
}

void Solver::vecsum(double *C, double a,const double *A, double b, const double *B,long int N)
{
	for(int I = 0; I < N; I++)
		C[I] = a*A[I] + b*B[I];
}

void Solver::vecscale(double *B, double a, double*A,long int N)
{
	for(int I = 0; I < N; I++)
		B[I] = a*A[I];
}

void Solver::matvecmul(double *y, const double *const *A, const double *x,long int N, int nNbr, const int *Nbr)
{
	for(int I = 0; I < N; I++)
	{
		double sum = 0.0;
		for(unsigned short t = 0; t < nNbr; t++)
			sum =  sum + A[t][I]*((I+Nbr[t] >= 0 && I + Nbr[t] < N)? x[I+Nbr[t]]:0.0);
		y[I] = sum;
	}
}

bool Solver::IsConverged(double *V, int N)
{
	double S = L1Norm(V,N);
	// cout<<S/(ResNorm + 1e-12)<<endl;
	if(S/(ResNorm+1e-12) < Tolerance)
		return true;
	return false;
}

void Solver::Initialize(long int *Map,long int l,long int u, long int iN)
{
	lNx = l; uNx = u; N = iN;
	Map21 = Map;

	Tolerance = TOL;
	alpha = 0.0;
	orf = 1.4;
	
	Allocate1Dvect(&Rn,N);

	if(Pcode == 1)
		Allocate2Dmat(&LU,5,N);

	if(Scode == 1 && Pcode == 0)
	{
		Allocate1Dvect(&phi_,N);
		Allocate1Dvect(&r_p,N);
		Allocate1Dvect(&p,N);
		Allocate1Dvect(&Ap,N);
		Allocate1Dvect(&r1,N);
		Allocate1Dvect(&s,N);
		Allocate1Dvect(&As,N);
		Allocate1Dvect(&s2,N);
		a1 = 1.0, a2 = 1.0, a3 = 1.0;
	}

	if(Scode == 1 && Pcode == 1)
	{
		Allocate1Dvect(&r_,N);
		Allocate1Dvect(&v,N);
		Allocate1Dvect(&w,N);
		Allocate1Dvect(&p,N);
		Allocate1Dvect(&p_,N);
		Allocate1Dvect(&s,N);
		Allocate1Dvect(&s_,N);
		rho1 = 1.0;
		rho2 = 1.0;
		Alpha = 1.0;
		Beta = 1.0;
		Omega = 1.0;
	}
}

Solver::~Solver()
{
	DeAllocate1Dvect(&Rn ,N);

	Map21 = NULL;

	if(LU) DeAllocate2Dmat(&LU,nNbr,N);

	if(Scode == 1 && Pcode == 0)
	{
		DeAllocate1Dvect(&phi_,N);
		DeAllocate1Dvect(&r_p,N);
		DeAllocate1Dvect(&p,N);
		DeAllocate1Dvect(&Ap,N);
		DeAllocate1Dvect(&r1,N);
		DeAllocate1Dvect(&s,N);
		DeAllocate1Dvect(&As,N);
		DeAllocate1Dvect(&s2,N);
	}

	if(Scode == 1 && Pcode == 1)
	{
		DeAllocate1Dvect(&r_,N);
		DeAllocate1Dvect(&v,N);
		DeAllocate1Dvect(&w,N);
		DeAllocate1Dvect(&p,N);
		DeAllocate1Dvect(&p_,N);
		DeAllocate1Dvect(&s,N);
		DeAllocate1Dvect(&s_,N);
	}
}

void Solver::ComputeResidual()
{
	matvecmul(Rn,A,phi,N,nNbr,Nbr);
	vecsum(Rn,1.0,b,-1.0,Rn,N);
}

int Solver::GSSORSolver()
{
	int iters = 0;
	long int P;
	while(1)
	{
		for(int I = lNx; I < uNx; I++)
		{
			for(P = Map21[I]+1; P < Map21[I+1]-1; P++)
			{
				double ai = 1.0, sum = 0.0;
				sum = A[0][P]*phi[P+Nbr[0]] + A[1][P]*phi[P+Nbr[1]] + A[2][P]*phi[P+Nbr[2]]
						+ A[3][P]*phi[P+Nbr[3]] + A[4][P]*phi[P+Nbr[4]];
				ai = A[2][P];
				phi[P] = phi[P] + orf*(b[P] - sum)/ai;
			}
		}
		ComputeResidual();
		if(iters == 0) ResNorm = L1Norm(Rn,N);
		iters++;
		cout<<Solver::L1Norm(Rn,N)<<endl;
		if(IsConverged(Rn,N))
			break;
	}
	// cout <<iters<<endl;
	return 0;
}

void Solver::SolveLU(double *x, double*b)
{
	long int I,J;
	unsigned short t;
	double sum = 0.0;
	//Forward Substitution
	for(I = lNx; I < uNx; I++)
		for(J = Map21[I]+1; J < Map21[I+1]-1; J++)
		{
			sum = 0.0;
			for(t = 0; t < nNbr/2; t++)
				sum += LU[t][J]*x[J+Nbr[t]];
			x[J] = (b[J] - sum)/LU[nNbr/2][J];
		}	

	//Backward Substitution
	for(I = uNx-1; I >= lNx; I--)
		for(J = Map21[I+1]-2; J>= Map21[I]+1; J--)
		{
			sum = 0.0;
			for(t = nNbr-1; t > nNbr/2; t--)
				sum += LU[t][J]*x[J+Nbr[t]];
			x[J] = (x[J] - sum);
		}
}

int Solver::PBiCGStabSolver()
{
	cout<<"Running preconditioned BiCGStabSolver"<<endl;

	if(Pcode == 1)
		ILU();

	ComputeResidual();
	copyvec(Rn, r_, N);
	ResNorm = L1Norm(Rn,N);
	int iters = 1;
	if((FinRes = L1Norm(Rn,N)/ResNorm) < 1e-5)
		goto retpoint0;

	if(fabs(ResNorm-0.0) < 1e-10)	ResNorm = 1.0;

	while(1)
	{
		rho1 = dotp(r_,Rn,N);
		if(fabs(rho1) < 1e-12)
		{
			FinRes = L1Norm(Rn,N)/ResNorm;
			goto retpoint1;
		}

		if(iters == 1)
			copyvec(Rn, p, N);
		else
		{
			Beta = (rho1/rho2)*(Alpha/Omega);
			vecsum(p,1.0,p,-Omega,v,N);
			vecsum(p,1.0,Rn,Beta,p,N);
		}

		SolveLU(p_,p);
		matvecmul(v,A,p_,N,nNbr,Nbr);
	
		Alpha = rho1/dotp(r_,v,N);
		vecsum(s,1.0,Rn,-Alpha,v,N);
		if((FinRes = L1Norm(s,N)/ResNorm) < 1e-5)
		{
			vecsum(phi,1.0,phi,-Alpha,p_,N);
			goto retpoint0;
		}

		SolveLU(s_,s);
		matvecmul(w,A,s_,N,nNbr,Nbr);
		
		Omega = dotp(w,s,N)/dotp(w,w,N);

		vecsum(phi,1.0,phi,Alpha,p_,N);
		vecsum(phi,1.0,phi,Omega,s_,N);
		

		vecsum(Rn,1.0,s,-Omega,w,N);

		rho2 = rho1;
		if((FinRes = L1Norm(Rn,N)/ResNorm) < 1e-5)
			goto retpoint0;

		if(abs(Omega) < 1e-12)
		{
			FinRes = L1Norm(Rn,N)/ResNorm;
			goto retpoint1;
		}

		iters++;
	}

	FinRes = L1Norm(Rn,N)/ResNorm;
	if(iters == inneriters)
		goto retpoint1;

	retpoint0:
	cout<<iters<<"	"<<FinRes<<endl;
	return 0;

	retpoint1:
	cout<<iters<<"	"<<FinRes<<endl;
	return 1;
}

void Solver::ILU()
{
	long int I,J;
	double A2[4];
	Allocate2Dmat(&LU,nNbr,N);
	for(I = lNx; I < uNx; I++)
	{
		for(J = Map21[I]+1; J < Map21[I+1]-1; J++)
		{
			LU[0][J] = A[0][J]/(1.0 + alpha*LU[3][J+Nbr[0]]);
			LU[1][J] = A[1][J]/(1.0 + alpha*LU[4][J+Nbr[1]]);
			
			A2[0] = alpha*LU[0][J]*LU[3][J+Nbr[0]];
			A2[1] = alpha*LU[1][J]*LU[4][J+Nbr[1]];
			A2[2] = LU[1][J]*LU[3][J+Nbr[1]];
			A2[3] = LU[0][J]*LU[4][J+Nbr[0]];

			LU[2][J] = A[2][J] + A2[0] + A2[1] + A2[2] + A2[3];
			LU[3][J] = (A[3][J]-A2[0])/LU[2][J];
			LU[4][J] = (A[4][J]-A2[1])/LU[2][J];
		}
	}
}

void Solver::Solve()
{
	if(Scode == 0)
		GSSORSolver();
	if(Scode == 1 && Pcode == 0)
		Solver::Getm() == 0? GSSORSolver() : BiCGStabSolver();
	if(Scode == 1 && Pcode == 1)
		PBiCGStabSolver();
}

int Solver::BiCGStabSolver()
{
	// cout<<"Running BiCGStabSolver"<<endl;
	ComputeResidual();
     copyvec(Rn, r1, N);
     copyvec(Rn, p, N);
     ResNorm = L1Norm(Rn,N);
     int iters = 0;
     while(iters < inneriters)
     {
		copyvec(Rn, r_p, N);
		matvecmul(Ap, A, p, N, nNbr, Nbr);
		a1 = dotp(r_p, r1, N) / dotp(Ap, r1, N);
		vecsum(s, 1.0, r_p, -a1, Ap, N);
		matvecmul(As, A, s, N, nNbr, Nbr);
		a2 = dotp(As, s, N) / dotp(As, As, N);
		vecsum(s2, a1, p, a2, s, N);
		vecsum(phi, 1.0, phi, 1.0, s2, N);
		vecsum(Rn, 1.0, s, -a2, As, N);
		if (IsConverged(Rn,N))
		   break;
		a3 = (a1 / a2) * dotp(Rn, r1, N) / dotp(r_p, r1, N);
		vecsum(s2, 1.0, p, -a2, Ap, N);
		vecsum(p, 1.0, Rn, a3, s2, N);
		if (fabs(dotp(Rn, r1, N)) < 1e-6)
		{
		   copyvec(Rn, r1, N);
		   copyvec(Rn, p, N);
		}
		//if (iters > 3000)
		//      break;
		iters++;
     }
     // cout<<iters<<"	"<<L1Norm(Rn,N)<<endl;
}