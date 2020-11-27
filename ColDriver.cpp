#include "ColDriver.h"

using namespace std;

unsigned short Driver::Time_Scheme=0;
unsigned short Driver::restart=0;
double Driver::dt=0.0;
double Driver::time=0.0;
unsigned long int Driver::ntime=0;
unsigned long int Driver::itime=0;


Driver::Driver()
{
	G = NULL;
	nDiag = 0, nvars = 0;
	A = NULL, S = NULL;
	phi = NULL, phi_ = NULL, phi__ = NULL;
	Ce = 0.0, Cp = 0.0, Cn = 0.0, D = 0.0;
	Lambda = NULL, Urf = NULL;
}

void Driver::Initialize(Geometry *inG, int in_nvars, int in_nDiag)
{
	if(inG == NULL)
	{
		cout<<"Geometry not initialized"<<endl;
		exit(0);
	}

	G = inG;
	nDiag = in_nDiag; nvars = in_nvars;
	Lambda = new double[nvars];
	Urf = new double[nvars];
	long int N = G->N;

	Allocate2Dmat(&A,nDiag,N);
	Allocate2Dmat(&S,nvars,N);
	Allocate2Dmat(&phi,nvars,N);

	
	NbrIdx = new int[5];
	//Initialized for 5 diagonals for now
	NbrIdx[0] = -1;
	NbrIdx[1] = -(G->CNy);
	NbrIdx[2] = 0;
	NbrIdx[3] = G->CNy;
	NbrIdx[4] = 1;

	if(Time_Scheme != 2)
		Allocate2Dmat(&phi_,nvars,N);
	if(Time_Scheme == 1)
		Allocate2Dmat(&phi__,nvars,N);	

	ntime = 0; time = 0.0;
}

Driver::~Driver()
{
	long int N = G->N;
	DeAllocate2Dmat(&A,nDiag,N);
	DeAllocate2Dmat(&S,nvars,N);
	DeAllocate2Dmat(&phi, nvars,N);

	if(phi_ != NULL) DeAllocate2Dmat(&phi_,nvars,N);
	if(phi__ != NULL) DeAllocate2Dmat(&phi__,nvars,N);

	ntime = 0; time = 0.0;
	nDiag = 0, nvars = 0;
	delete[] Lambda, Urf;
	Lambda = NULL, Urf = NULL;
	G = NULL;
	cout<<"Driver destroyed"<<endl;
}

void Driver::StoreMat(string fn)
{
	fstream f;
	f.open(fn,ios::out);
	long int P; long int CNx = G->CNx, CNy = G->CNy;

	for(long int I = 1; I < CNx-1; I++)
	{	
		P = G->Map21[I];
		for(long int J = 1; J < CNy-1; J++)
			f<<I<<"	"<<J<<"	"<<A[0][P+J]<<"	"<<A[1][P+J]<<"	"<<A[2][P+J]<<"	"<<A[3][P+J]<<"	"<<A[4][P+J]<<endl;
	}
	f.close();
}

void Driver::StoreS(string fn, int varc)
{
	fstream f;
	f.open(fn,ios::out);
	long int P; long int CNx = G->CNx, CNy = G->CNy;

	for(long int I = 1; I < CNx-1; I++)
	{	
		P = G->Map21[I];
		for(long int J = 1; J < CNy-1; J++)
			f<<I<<"	"<<J<<"	"<<S[varc][P+J]<<endl;
	}
	f.close();

}

MomDriver::MomDriver(double irho, double imu) : rho(irho), mu(imu), bodyforcing(0)
{
	vF = NULL;
	p = NULL, mdotE = NULL, mdotN = NULL;
	Ap = NULL, pgrad = NULL;	
}

MomDriver::MomDriver(double irho, double imu, double igx, double igy, double ibeta) : rho(irho), mu(imu),
							gx(igx), gy(igy), beta(ibeta), bodyforcing(1)
{
	vF = NULL;
	p = NULL, mdotE = NULL, mdotN = NULL;
	Ap = NULL, pgrad = NULL;
}

void MomDriver::Initialize(Geometry *inG, int nDim, int nDiag)
{
	Time_Scheme = 2;
	inners = 80;
	Tol = 0.02;
	Driver::Initialize(inG, nDim, nDiag);
	Allocate2Dmat(&Ap,nDim,G->N);
	Allocate2Dmat(&pgrad,nDim,G->N);
	Sol = new Solver();
	Sol->SetScode(1);
	Sol->SetPcode(0);
	Sol->Initialize(G->Map21, 1, G->CNx-1, G->N);
	Sol->SetA(A, NbrIdx, 5);
	Sol->SetnInner(inners);
	Sol->SetTol(Tol);
	// cout<<Ap[U]<<"	"<<Ap[V]<<endl;
	double Ur = max(Solver::LinfNorm(phi[U],G->N),Solver::LinfNorm(phi[V],G->N));
	double Lr = min(G->x[G->Nx-1]-G->x[0],G->y[G->Nx-1]-G->y[0]);
	cout<<"Reynolds Number = "<<rho*Ur*Lr/mu<<endl;
	if(bodyforcing == 1)
		cout<<"Rayleigh Number = "<<rho*rho*max(gx,gy)*beta*(1.0)*Lr*Lr*Lr*0.1/(mu*mu)<<endl;

	Lambda[U] = 1.0, Urf[U] = 0.8;
	Lambda[V] = 1.0, Urf[V] = 0.8;

	cout<<"Momentum Driver Setup"<<endl;
}

void MomDriver::Solve()
{
	if(Driver::GetTimeScheme() !=2)
	{
		if(Driver::GetTimeScheme() == 1)
		{
			Solver::copyvec(phi_[U],phi__[U],G->N);
			Solver::copyvec(phi_[V],phi__[V],G->N);	
		}
		Solver::copyvec(phi[U],phi_[U],G->N);
		Solver::copyvec(phi[V],phi_[V],G->N);	
	}

	ComputeConDiff();
	ComputePGrad();
	if(bodyforcing)
		ComputeSource();
	ComputeTime();
	SetBC();
	
	SetEquation(U);
	Sol->Setb(S[U]);
	Sol->Setphi(phi[U]);
	Sol->Solve();
	MomRes[U] = Sol->GetRes();
	SetEquation(V);
	Sol->Setb(S[V]);
	Sol->Setphi(phi[V]);
	Sol->Solve();
	MomRes[V] = Sol->GetRes();
}

void MomDriver::ComputeConDiff()
{
	long int Nc = G->N;
	for(long int P = 0; P < Nc; P++)
	{
		for(unsigned short i = 0; i < nvars; i++)
		{
			Ap[i][P] = 0.0;
			S[i][P] = 0.0;
		}
	}

	long int CNx = G->CNx, CNy = G->CNy;
	double gxE, gxP, gyN, gyP;
	double dxE, dyN, Se, Sn;
	long int P, E, N;
	double UDS_[2],CDS_[2];


	//Looping through all the internal east faces of control volumes
	for(long int I = 1; I <  CNx-2; I++)
	{
		gxE = G->gx[I]; gxP = 1.0 - gxE;
		dxE = G->X[I+1] - G->X[I];

		P = G->Map21[I]; E = P + CNy;
		for(long int J = 1; J < CNx-1; J++)
		{
			P = P+1;
			E = E+1;

			Se = G->y[J] - G->y[J-1];
			D = mu*Se/dxE;

			Ce = min(mdotE[P],0.0); Cp = max(mdotE[P],0.0);

			A[3][P] = Ce-D;
			A[1][E] = -Cp-D; 

			//Deferred Correction sources
			UDS_[U] = Ce*phi[U][E] + Cp*phi[U][P];
			UDS_[V] = Ce*phi[V][E] + Cp*phi[V][P]; 

			CDS_[U] = mdotE[P]*(phi[U][P]*gxP + phi[U][E]*gxE);
			CDS_[V] = mdotE[P]*(phi[V][P]*gxP + phi[V][E]*gxE);

			S[U][P] += Lambda[U]*(-CDS_[U] + UDS_[U]);
			S[U][E] -= Lambda[U]*(-CDS_[U] + UDS_[U]);

			S[V][P] += Lambda[V]*(-CDS_[V] + UDS_[V]);
			S[V][E] -= Lambda[V]*(-CDS_[V] + UDS_[V]);
		}
	}

	//Looping through all the internal north faces of control volumes
	for(long int J = 1; J < CNy-2; J++)
	{
		gyN = G->gy[J]; gyP = 1.0-gyN;
		dyN = G->Y[J+1]-G->Y[J];

		for(long int I = 1; I < CNx-1; I++)
		{
			P = G->Map21[I]+J; N = P+1;
			Sn = G->x[I]-G->x[I-1];

			D = mu*Sn/dyN;

			Cn = min(mdotN[P],0.0); Cp = max(mdotN[P],0.0);

			A[0][N] = -Cp-D;
			A[4][P] = Cn-D;

			//Deferred correction
			UDS_[U] = Cn*phi[U][N] + Cp*phi[U][P];
			UDS_[V] = Cn*phi[V][N] + Cp*phi[V][P];

			CDS_[U] = mdotN[P]*(gyP*phi[U][P] + gyN*phi[U][N]);
			CDS_[V] = mdotN[P]*(gyP*phi[V][P] + gyN*phi[V][N]);

			S[U][P] += Lambda[U]*(-CDS_[U] + UDS_[U]);
			S[U][N] -= Lambda[U]*(-CDS_[U] + UDS_[U]);
			S[V][P] += Lambda[V]*(-CDS_[V] + UDS_[V]);
			S[V][N] -= Lambda[U]*(-CDS_[V] + UDS_[V]);
		}
	}
}

void MomDriver::ComputePGrad()
{
	double Dx, Dy;
	long int P, e,n,w,s;
	long int CNx = G->CNx, CNy = G->CNy;
	double gxe, gyn, gxw, gys;
	double pe, pw, pn, ps;

	for(long int I = 1; I < CNx-1; I++)
	{
		Dx = G->x[I]-G->x[I-1];
		gxe = G->gx[I];gxw = G->gx[I-1];
		P = G->Map21[I];
		for(long int J = 1; J < CNy-1; J++)
		{
			Dy = G->y[J]-G->y[J-1];
			P = P+1;
			e = P + CNy; w = P - CNy;
			n = P + 1; s = P - 1;
			gyn = G->gy[J]; gys = G->gy[J-1];

			pe = gxe*p[e] + (1.0-gxe)*p[P];
			pw = gxw*p[P] + (1.0-gxw)*p[w];
			pn = gyn*p[n] + (1.0-gyn)*p[P];
			ps = gys*p[P] + (1.0-gys)*p[s];

			pgrad[U][P] = (pe-pw)*Dy;
			pgrad[V][P] = (pn-ps)*Dx;

			S[U][P] += (-pgrad[U][P]);
			S[V][P] += (-pgrad[V][P]);
		}
	}
}

void MomDriver::ComputeTime()
{
	if(Time_Scheme == 2)
		return;
		

	long int P;
	long int CNx = G->CNx, CNy = G->CNy;
	double Dx, Dy, Vc, At;
	double dt = Driver::GetTimeStep();
	for(int I = 1; I < CNx-1; I++)
	{
		P = G->Map21[I];
		Dx = G->x[I]-G->x[I-1];
		for(int J = 1; J < CNy-1; J++)
		{
			P = P+1;
			Dy = G->y[J]-G->y[J-1];
			Vc = Dx*Dy;
			At = rho*Vc/dt;
			if(Time_Scheme == 0)
			{
				S[U][P] += At*phi_[U][P];
				S[V][P] += At*phi_[V][P];
				Ap[U][P] += At;Ap[V][P] += At;
			}
			if(Time_Scheme == 1)
			{
				S[U][P] += At*0.5*(4*phi_[U][P]-phi__[U][P]);
				S[V][P] += At*0.5*(4*phi_[V][P]-phi__[V][P]);
				Ap[U][P] += At*1.5; Ap[V][P] += At*1.5;
			}
		}
	}
}

void MomDriver::SetEquation(int varc)
{
	long int CNx = G->CNx, CNy = G->CNy;
	fstream f;
	char fn[30];
	sprintf(fn,"Ar_%d.dat",varc);
	// f.open(fn,ios::out);
	for(long int I = 1; I < CNx-1; I++)
	{
		long int P = G->Map21[I];
		for(long int J = 1; J < CNy-1; J++)
		{
			P = P+1;
			A[2][P] = (-(A[0][P] + A[1][P] + A[3][P] + A[4][P]) + Ap[varc][P])/Urf[varc];
			S[varc][P] += A[2][P]*phi[varc][P]*(1.0 - Urf[varc]);
			Ap[varc][P] = 1.0/A[2][P];
			// f<<I<<"	"<<J<<"	"<<Ap[varc][P]<<"	"<<A[2][P]<<endl;
		}
	}
	// f.close();
}

void MomDriver::ComputeSource()
{
	long int P;
	int CNx = G->CNx, CNy = G->CNy;
	double Dx, Dy, Vc;
	for(long int I = 1; I < CNx-1; I++)
	{
		P = G->Map21[I];
		Dx = G->x[I] - G->x[I-1];
		for(int J = 1; J < CNy-1; J++)
		{
			P = P+1;
			Dy = G->y[J] - G->y[J-1];
			Vc = Dx*Dy;
			S[U][P] += -rho*gx*beta*Vc*(vF[P]-1.0);
			S[V][P] += -rho*gy*beta*Vc*(vF[P]-1.0);
		}
	}
}

void MomDriver::SetBC()
{
	SetMomBC(U);
	SetMomBC(V);
}

void MomDriver::SetMomBC(int varc)
{
	int CNx = G->CNx, CNy = G->CNy;
	long int P;
	
	//South Boundary
	double n[2] = {0.0,-1.0};
	for(long int I = 1; I < CNx-1; I++)
	{
		P = G->Map21[I]+1;
		D = mu*(G->x[I] - G->x[I-1])*(1-n[varc]*n[varc])/((G->Y[0] - G->Y[1])*n[1]);
		Ap[varc][P] += D;
		S[varc][P] += (D)*phi[varc][P-1];
	}


	//West Boundary
	n[0] = -1.0, n[1] = 0.0;
	P = G->Map21[1];
	for(long int J = 1; J < CNy-1; J++)
	{
		P = P+1;
		D = mu*(G->y[J] - G->y[J-1])*(1-n[varc]*n[varc])/((G->X[0]-G->X[1])*n[0]);
		Ap[varc][P] += D;
		// if(varc == V)
		// 	cout<<P<<"	"<<Ap[varc][P]<<"	"<<D<<endl;
		S[varc][P] += (D)*phi[varc][P-CNy];
	}

	//North Boundary
	n[0] = 0.0, n[1]= 1.0;
	for(long int I = 1; I < CNx-1; I++)
	{
		P = G->Map21[I] + CNy-2;
		D = mu*(G->x[I] - G->x[I-1])*(1-n[varc]*n[varc])/((G->Y[CNy-1] - G->Y[CNy-2])*n[1]);
		Ap[varc][P] += D;
		S[varc][P] += (D)*phi[varc][P+1];
	}

	//East Boundary
	n[0] = 1.0, n[1] = 0.0;
	P = G->Map21[CNx-2];
	for(long int J = 1; J < CNy-1; J++)
	{
		P = P + 1;
		D = mu*(G->y[J] - G->y[J-1])*(1-n[varc]*n[varc])/((G->X[CNx-1]-G->X[CNx-2])*n[0]);
		Ap[varc][P] += D;
		// if(varc == V)
		// 	cout<<P<<"	"<<Ap[varc][P]<<"	"<<D<<endl;
		S[varc][P] += (D)*phi[varc][P+CNy];
	}
}

void MomDriver::Store()
{
	int CNx = G->CNx, CNy = G->CNy;
	long int P;
	fstream fu;
	char fn[50], gn[50], hn[50];
	
	sprintf(fn,"Velocity_%ld.dat",Driver::GetiTime());
	sprintf(gn,"Streamfunction_%ld",Driver::GetiTime());
	

	fu.open(fn,ios::out);
	double dl = pow(G->X[1]-G->X[0],2) + pow(G->Y[1]-G->Y[0],2);
	dl = sqrt(dl);
	double u,v, dx, dy, theta;
	for(long int I = 0; I < CNx; I++)
	{
		P = G->Map21[I];
		for(long int J = 0; J < CNy; J++)
		{
			u = phi[U][P+J]; v = phi[V][P+J];
			if(fabs(u) < TOL && fabs(v) < TOL){dx = 0.0; dy = 0.0;}
			else{dx = dl*u/sqrt(u*u+v*v); dy = dl*v/sqrt(u*u+v*v);}
			fu<<G->X[I]<<"	"<<G->Y[J]<<"	"<<u<<"	"<<v<<"	"<<dx<<"	"<<dy<<endl;
		}
		fu<<endl;
	}
	fu.close();

	double **sf;
	Allocate2Dmat(&sf, CNx, CNy);
	fu.open(gn,ios::out);
	long int I, J;
	for(J = 1; J < CNy-1; J++)
		sf[0][J] = sf[0][J-1] + mdotE[G->Map21[0]+J];

	for(I = 1; I < CNx; I++)
	{
		P = G->Map21[I];
		sf[I][0] = sf[I-1][0] - mdotN[P];
		for(int J = 1; J < CNy; J++)
			sf[I][J] = sf[I][J-1] + mdotE[P+J];
	}

	double sfmax = 0.0;
	double sfmin = 0.0;
	for(I = 0; I < CNx; I++)
	{
		for(int J = 0; J < CNy; J++)
		{
			sfmax = max(sf[I][J],sfmax);
			sfmin = min(sf[I][J],sfmin);
			fu<<G->X[I]<<"	"<<G->Y[J]<<"	"<<sf[I][J]<<endl;
		}
		fu<<endl;
	}
	fu.close();

	if(Driver::GetTimeScheme() == 2)
	{
		fu.open("Sf_convergence.dat",ios::app);
		fu<<(G->Nx-1)*(G->Ny-1)<<"	"<<sfmax<<"	"<<sfmin<<endl;
		fu.close();
	}
		
	cout<<"Max Sf value = "<<sfmax<<endl;
	cout<<"Min Sf value = "<<sfmin<<endl;
	DeAllocate2Dmat(&sf,CNx,CNy);



}

PresDriver::PresDriver(double irho): rho(irho)
{
	p = NULL, mdotE = NULL, mdotN = NULL;
	pgradx = NULL, pgrady = NULL;
}

void PresDriver::Initialize(Geometry *inG, int nvar, int nDiag)
{
	Driver::Initialize(inG,nvar,nDiag);
	Allocate1Dvect(&p,G->N);
	Allocate1Dvect(&mdotE,G->N);
	Allocate1Dvect(&mdotN,G->N);
	Sol = new Solver();
	Sol->SetScode(1);
	Sol->SetPcode(0);
	Sol->Initialize(G->Map21, 1, G->CNx-1, G->N);
	Sol->SetA(A, NbrIdx, 5);
	Sol->SetTol(0.01);
	Sol->SetnInner(900);
	Lambda[0] = 0.0; Urf[0] = 0.3;
	cout<<"Pressure Driver Setup"<<endl;
}

void PresDriver::Init()
{
	if(Driver::Getrestart() == 1)
	{
		long int P = 0;
		fstream f,g,h;
		f.open("U.dat",ios::in);g.open("V.dat",ios::in);h.open("P.dat",ios::in);
		for(long int I = 1; I < G->CNx-1; I++)
		{
			P = G->Map21[I];
			for(long int J = 1; J < G->CNy-1; J++)
			{
				f>>uf[P+J]; g>>vf[P+J]; h>>p[P+J];
			}

		}
		f.close(); g.close(); h.close();

		double gxE, gyN, Se, Sn;
		for(long int I = 1; I < G->CNx-1; I++)
		{
			P = G->Map21[I];
			gxE = G->gx[I];
			Sn = G->x[I]-G->x[I-1];
			for(long int J = 1; J < G->CNy-1; J++)
			{
				gyN = G->gy[J];
				Se = G->y[J]-G->y[J-1];
				mdotE[P+J] = rho*Se*(gxE*uf[P+G->CNy]+(1.0-gxE)*uf[P]);
				mdotN[P+J] = rho*Sn*(gyN*vf[P+G->CNy]+(1.0-gyN)*vf[P]);
			}
		}
	}
		

	SetBC();
}

PresDriver::~PresDriver()
{
	DeAllocate1Dvect(&mdotN,G->N);
	DeAllocate1Dvect(&mdotE,G->N);
	DeAllocate1Dvect(&p,G->N);
}

void PresDriver::Solve()
{
	ComputeConDiff();
	SetBC();
	SetEquation();
	// StoreMat("Ap.dat");
	// StoreS("bp.dat",0);
	Sol->Setb(S[0]);
	// cout<<"Residual Sum = "<<Solver::sumvec(S[0],G->N)<<endl;
	Sol->Setphi(phi[0]);
	Sol->Solve();
	UpdateBoundaryPres(0);
	Corrections();
	// cout<<"Max u Vel = "<<Solver::LinfNorm(uf,G->N)<<endl;
	// cout<<"Max v Vel = "<<Solver::LinfNorm(vf,G->N)<<endl;
	Pres = Sol->GetRes();

	//Mom[U]->Solve(A,S[U]);
}

void PresDriver::ComputeConDiff()
{
	//A pressure weighted Rhie-Chow interpolation scheme is used for pressure velocity coupling
	

	//Loop over internal east faces of CV's
	long int CNx = G->CNx, CNy = G->CNy;
	double dxE, dyN, gxE, gyN;
	double Vnf, dp, Af, bf, Sf;
	long int P, E, N;

	for(long int I = 1; I < CNx-2; I++)
	{
		dxE = G->X[I+1] - G->X[I];
		gxE = G->gx[I];

		P = G->Map21[I]; E = P + CNy;
		for(long int J = 1; J < CNy-1; J++)
		{
			P = P+1; E = E + 1;
			Sf = G->y[J]-G->y[J-1];

			dp = -(p[E]-p[P])*Sf;
			Af = gxE*Apu[E] + (1.0-gxE)*Apu[P];
			bf = gxE*(u[E]/Apu[E] + pgradx[E]) + (1.0-gxE)*(u[P]/Apu[P] + pgradx[P]);
			Vnf = (bf + dp)*Af;

			D = rho*Sf;
			mdotE[P] = D*Vnf;
			A[3][P] = -D*Af*Sf;
			A[1][E] = A[3][P];
		}
	}


	//Loop through internal North faces of CVs
	for(long int I = 1; I < CNx-1; I++)
	{
		Sf = G->x[I] - G->x[I-1];
		P = G->Map21[I]; N = P+1;
		for(long int J = 1; J < CNy-2; J++)
		{
			P = P+1; N = N+1;
			dyN = G->Y[J+1]-G->Y[J];
			gyN = G->gy[J];

			dp = -(p[N]-p[P])*Sf;
			Af = gyN*Apv[N] + (1.0-gyN)*Apv[P];
			bf = gyN*(v[N]/Apv[N] + pgrady[N]) + (1.0-gyN)*(v[P]/Apv[P] + pgrady[P]);
			Vnf = (bf + dp)*Af;

			D = rho*Sf;
			mdotN[P] = D*Vnf;
			A[4][P] = -D*Af*Sf;
			A[0][N] = A[4][P];
		}
	}
}

void PresDriver::SetBC()
{
	SetPresBC();

	//Velocity BC for Lid Driven Cavity
	for(int I = 1; I < G->CNx-1; I++)
	{
		uf[G->Map21[I] + G->CNy-1] = 1.0;
		mdotE[G->Map21[I] + G->CNy-1] = rho*uf[G->Map21[I] + G->CNy-1]*(G->Y[G-> CNy-1] - G->y[G->Ny-1]);
	}
}

void PresDriver::SetPresBC()
{
	long int P;

	//South Boundary & North Boundary
	for(long int I = 1; I < G->CNx-1; I++)
	{
		mdotN[G->Map21[I]] = 0.0;
		mdotN[G->Map21[I]+G->CNy-1] = 0.0;
	}

	//West Boundary & East Boundary
	for(long int J = 1; J < G->CNy-1; J++)
	{
		mdotE[J] = 0.0;
		P = G->Map21[G->CNx-1];
		mdotE[P+J] = 0.0;
	}	

	Mdot = 0.0;
}

void PresDriver::SetEquation(int varc)
{
	long int P;
	int CNx = G->CNx, CNy = G->CNy;
	fstream f;
	for(long int I = 1; I < CNx-1; I++)
	{
		P = G->Map21[I];
		for(long int J = 1; J < CNy-1; J++)
		{
			P = P+1;
			S[varc][P] = -(mdotE[P] - mdotE[P-CNy] + mdotN[P]-mdotN[P-1]);
			A[2][P] = -(A[0][P] + A[1][P] + A[3][P] + A[4][P]);
			MdotC += S[varc][P];
		}
	}
}


void PresDriver::UpdateBoundaryPres(int varc)
{
	long int P;
	double g;
	double *lp;
	if(varc == 0)	lp = phi[0];
	if(varc == 1)	lp = p;

	long int CNx = G->CNx, CNy = G->CNy;

	//South Boundary & North Boundary
	for(long int I = 1; I < CNx-1; I++)
	{
		P = G->Map21[I];
		g = G->gy[1];
		lp[P] = lp[P+1] + g*(lp[P+1]-lp[P+2]);

		P = G->Map21[I] + CNy-1;
		g = G->gy[CNy-3];
		lp[P] = lp[P-1] + (1.0-g)*(lp[P-1]-lp[P-2]);
	}

	//West Boundary & East Boundary
	for(long int J = 1; J < CNy-1; J++)
	{
		P = J;
		g = G->gx[1];
		lp[P] = lp[P+CNy] + g*(lp[P+CNy]-lp[P+2*CNy]);

		P = G->Map21[CNx-1]+J;
		g = G->gx[CNx-3];
		lp[P] = lp[P-CNy] + (1.0-g)*(lp[P-CNy]-lp[P-2*CNy]);
	}

	lp = NULL;
}

void PresDriver::Corrections()
{
	long int P;
	
	//Correct Interal cell face mass fluxes
	long int CNx = G->CNx, CNy = G->CNy;
	for(long int I = 1; I < CNx-2; I++)
	{
		P = G->Map21[I];
		for(long int J = 1; J < CNy-1; J++)
			mdotE[P+J] = mdotE[P+J]+A[3][P+J]*(phi[0][P+J+CNy] - phi[0][P+J]);
	}

	for(long int I = 1; I < CNx-1; I++)
	{
		P = G->Map21[I];
		for(int J = 1; J < CNy-2; J++)
			mdotN[P+J] = mdotN[P+J]+A[4][P+J]*(phi[0][P+J+1] - phi[0][P+J]);
	}

	fstream f;
	Mdot = 0.0;
	// f.open("Mass_Check.dat",ios::out);
	for(int I = 1; I < CNx-1; I++)
	{
		P = G->Map21[I];
		for(long int J = 1; J < CNy-1; J++)
		{
			Mdot +=	(mdotE[P+J] + mdotN[P+J] - mdotE[P+J-CNy] - mdotN[P+J-1]);
			// f<<mdotE[P+J]<<"	"<<mdotN[P+J]<<"	"<<mdotE[P+J-CNy]<<"	"<<mdotN[P+J-1]<<"	"<<Mdot<<endl;
		}
	}
	// f.close();

	double pce, pcw, pcn, pcs;
	double gxE, gxP, gyN, gyP;
	double Dx, Dy;
	//Correct cell Center velocities and pressures
	pcref = phi[0][G->Map21[1] + 1];
	for(long int I = 1; I < CNx-1; I++)
	{
		P = G->Map21[I];
		gxE = G->gx[I]; gxP = G->gx[I-1];
		Dx = G->x[I]- G->x[I-1];
		for(long int J = 1; J < CNy-1; J++)
		{
			gyN = G->gy[J]; gyP = G->gy[J-1];
			Dy = G->y[J] - G->y[J-1];

			p[P+J] = p[P+J] + Urf[0]*(phi[0][P+J]-pcref);

			pce = (1.0-gxE)*phi[0][P+J] + gxE*phi[0][P+CNy+J];
			pcw = (1.0-gxP)*phi[0][P-CNy+J] + gxP*phi[0][P+J];
			pcn = (1.0-gyN)*phi[0][P+J] + gyN*phi[0][P+1+J];
			pcs = (1.0-gyP)*phi[0][P+J-1] + gyP*phi[0][P+J];

			uf[P+J] = uf[P+J] - (pce-pcw)*Dy*Apu[P+J];
			vf[P+J] = vf[P+J] - (pcn-pcs)*Dx*Apv[P+J];
		}
	}

	//Update Boundary Pressures
	UpdateBoundaryPres(1);
}

void PresDriver::Store()
{
	long int CNx = G->CNx, CNy = G->CNy;
	long int P;
	fstream fu;
	char fn[50];
	sprintf(fn,"Pressure_%ld.dat",Driver::GetiTime());
	
	fu.open(fn,ios::out);
	for(long int I = 0; I < CNx; I++)
	{
		P = G->Map21[I];
		for(long int J = 0; J < CNy; J++)
			fu<<G->X[I]<<"	"<<G->Y[J]<<"	"<<p[P+J]<<endl;
		fu<<endl;
	}
	fu.close();
}

ScalDriver::ScalDriver(double irho, double imu, double iPr, double iCp) : rho(irho), mu(imu), Pr(iPr), Cp(iCp)
{
	mdotE = NULL, mdotN = NULL;
}

ScalDriver::~ScalDriver()
{
	mdotE = NULL, mdotN = NULL;
	rho = 0, mu = 0, Pr = 0;
}

void ScalDriver::Initialize(Geometry *inG, int nDim, int nDiag)
{
	Time_Scheme = 2;
	Driver::Initialize(inG, nDim, nDiag);
	Sol = new Solver();
	Sol->SetScode(1);
	Sol->SetPcode(0);
	Sol->Initialize(G->Map21, 1, G->CNx-1, G->N);
	Sol->SetA(A, NbrIdx, 5);
	Sol->SetnInner(50);
	Sol->SetTol(0.05);
	// cout<<Ap[U]<<"	"<<Ap[V]<<endl;
	Lambda[T] = 1.0, Urf[T] = 0.9;

	cout<<"Scalar Driver Setup"<<endl;
}

void ScalDriver::Init()
{
	if(Driver::Getrestart() == 1)
	{
		long int P = 0;
		fstream f;
		f.open("T.dat",ios::in);
		for(long int I = 1; I < G->CNx-1; I++)
		{
			P = G->Map21[I];
			for(long int J = 1; J < G->CNy-1; J++)
				f>>phi[0][P+J];
		}
		f.close();
	}
	SetBC();
}

void ScalDriver::ComputeConDiff()
{
	long int Nc = G->N;
	long int CNx = G->CNx, CNy = G->CNy;
	double gxE, gxP, gyN, gyP;
	double dxE, dyN, Se, Sn;
	long int P, E, N;
	double UDS_[2],CDS_[2];


	//Looping through all the internal east faces of control volumes
	for(long int I = 1; I <  CNx-2; I++)
	{
		gxE = G->gx[I]; gxP = 1.0 - gxE;
		dxE = G->X[I+1] - G->X[I];

		P = G->Map21[I]; E = P + CNy;
		for(long int J = 1; J < CNx-1; J++)
		{
			P = P+1;
			E = E+1;

			Se = G->y[J] - G->y[J-1];
			D = (mu/Pr)*Se/dxE;

			Ce = min(mdotE[P],0.0); Cp = max(mdotE[P],0.0);

			A[3][P] = Ce-D;
			A[1][E] = -Cp-D; 

			//Deferred Correction sources
			UDS_[T] = Ce*phi[T][E] + Cp*phi[T][P];
			CDS_[T] = mdotE[P]*(phi[T][P]*gxP + phi[T][E]*gxE);
	
			S[T][P] += Lambda[T]*(-CDS_[T] + UDS_[T]);
			S[T][E] -= Lambda[T]*(-CDS_[T] + UDS_[T]);
		}
	}

	//Looping through all the internal north faces of control volumes
	for(long int J = 1; J < CNy-2; J++)
	{
		gyN = G->gy[J]; gyP = 1.0-gyN;
		dyN = G->Y[J+1]-G->Y[J];

		for(long int I = 1; I < CNx-1; I++)
		{
			P = G->Map21[I]+J; N = P+1;
			Sn = G->x[I]-G->x[I-1];

			D = (mu/Pr)*Sn/dyN;

			Cn = min(mdotN[P],0.0); Cp = max(mdotN[P],0.0);

			A[0][N] = -Cp-D;
			A[4][P] = Cn-D;

			//Deferred correction
			UDS_[T] = Cn*phi[T][N] + Cp*phi[T][P];

			CDS_[T] = mdotN[P]*(gyP*phi[T][P] + gyN*phi[T][N]);

			S[T][P] += Lambda[T]*(-CDS_[T] + UDS_[T]);
			S[T][N] -= Lambda[T]*(-CDS_[T] + UDS_[T]);
		}
	}
}

void ScalDriver::SetBC()
{
	long int Pw = G->Map21[0], Pe = G->Map21[G->CNx-1];
	if(Solver::Getm() == 0)
		for(long int J = 0; J < G->CNy; J++)
		{
			phi[T][Pw + J] = 0.0;
			phi[T][Pe + J] = 1.0;
		}
	SetScalBC();
}

void ScalDriver::SetScalBC()
{
	int CNx = G->CNx, CNy = G->CNy;
	long int P1,P2;
	double g = 0.0;

	
	//South Boundary & North Boundary - Adiabatic walls
	for(long int I = 1; I < CNx-1; I++)
	{
		P1 = G->Map21[I];
		g = G->gy[1];
		phi[T][P1] = phi[T][P1+1] + g*(phi[T][P1+1]-phi[T][P1+2]);

		P1 = P1 + CNy-1;
		g = G->gy[CNy-3];
		phi[T][P1] = phi[T][P1-1] + (1.0-g)*(phi[T][P1-1]-phi[T][P1-2]);
	}


	//West Boundary & East Boundary - Isothermal walls
	P1 = G->Map21[1]; P2 = G->Map21[CNx-2];
	for(long int J = 1; J < CNy-1; J++)
	{
		D = (mu/Pr)*(G->y[J] - G->y[J-1])/(-G->X[0]+G->X[1]);
		A[2][P1+J] += D;
		S[T][P1+J] += (D)*phi[T][P1+J-CNy];

		D = (mu/Pr)*(G->y[J] - G->y[J-1])/(G->X[CNx-1]-G->X[CNx-2]);
		A[2][P2+J] += D;
		S[T][P2+J] += (D)*phi[T][P2+J+CNy];
	}
}

void ScalDriver::ComputeTime()
{
	if(Time_Scheme == 2)
		return;
		

	long int P;
	long int CNx = G->CNx, CNy = G->CNy;
	double Dx, Dy, Vc, At;
	double dt = Driver::GetTimeStep();
	for(int I = 1; I < CNx-1; I++)
	{
		P = G->Map21[I];
		Dx = G->x[I]-G->x[I-1];
		for(int J = 1; J < CNy-1; J++)
		{
			P = P+1;
			Dy = G->y[J]-G->y[J-1];
			Vc = Dx*Dy;
			At = rho*Vc/dt;
			if(Time_Scheme == 0)
			{
				S[T][P] += At*phi_[T][P];
				A[2][P] += At;
			}
			if(Time_Scheme == 1)
			{
				S[T][P] += At*0.5*(4*phi_[T][P]-phi__[T][P]);
				A[2][P] += At*1.5;
			}
		}
	}
}

void ScalDriver::SetEquation(int varc)
{
	long int CNx = G->CNx, CNy = G->CNy;
	// fstream f;
	// char fn[30];
	// sprintf(fn,"Ar_%d.dat",varc);
	for(long int I = 1; I < CNx-1; I++)
	{
		long int P = G->Map21[I];
		for(long int J = 1; J < CNy-1; J++)
		{
			P = P+1;
			A[2][P] = (-(A[0][P] + A[1][P] + A[3][P] + A[4][P]) + A[2][P])/Urf[T];
			S[T][P] += A[2][P]*phi[T][P]*(1.0 - Urf[T]);
		}
	}
}

void ScalDriver::Solve()
{
	ComputeConDiff();
	SetBC();
	SetEquation();
	// StoreMat("AT.dat");
	// StoreS("bT.dat",0);
	Sol->Setb(S[T]);
	// cout<<"Residual Sum = "<<Solver::sumvec(S[T],G->N)<<endl;
	Sol->Setphi(phi[T]);
	Sol->Solve();
	Reset1Dvect(&A[2],G->N);
	Reset1Dvect(&S[T],G->N);
	// cout<<"Max T = "<<Solver::LinfNorm(phi[T],G->N)<<endl;
	Temp = Sol->GetRes();
}

void ScalDriver::Store()
{
	long int CNx = G->CNx, CNy = G->CNy;
	long int P;
	fstream fu;
	char fn[50], gn[50], hn[50];
	sprintf(fn,"Temperature_%ld.dat",Driver::GetiTime());
	
	fu.open(fn,ios::out);
	for(long int I = 0; I < CNx; I++)
	{
		P = G->Map21[I];
		for(long int J = 0; J < CNy; J++)
			fu<<G->X[I]<<"	"<<G->Y[J]<<"	"<<phi[T][P+J]<<endl;
		fu<<endl;
	}
	fu.close();

	fu.open("Nu_Walls.dat",ios::out);
	double Nu_L =  0.0, Nu_R = 0.0, f_L = 0.0, f_R = 0.0, Net = 0.0;
	for(long int I = 0; I < CNy; I++)
	{
		f_L = (mu/Pr)*(phi[T][G->Map21[1]+I]-phi[T][I])/(G->X[1]-G->X[0]);
		f_R = -(mu/Pr)*(phi[T][G->Map21[CNx-1]+I]-phi[T][G->Map21[CNx-2]+I])/(G->X[CNx-1]-G->X[CNx-2]);
		Nu_L += f_L* (I < CNy-2 ? (G->y[I+1]-G->y[I]):0.0);
		Nu_R += f_R* (I < CNy-2 ? (G->y[I+1]-G->y[I]):0.0);
		fu<<G->Y[I]<<"	"<<f_L<<"	"<<f_R<<endl;
		Net += f_L + f_R;
	}
	fu.close();

	Nu_L = Nu_L/(G->y[G->Ny-1] - G->y[0]);
	Nu_R = Nu_R/(G->y[G->Ny-1] - G->y[0]);
	fu.open("Nu_avg.dat",ios::app);
	fu<<(G->Nx-1)*(G->Ny-1)<<"	"<<Nu_L<<"	"<<Nu_R<<"	"<<Net<<endl;
	fu.close();

}





