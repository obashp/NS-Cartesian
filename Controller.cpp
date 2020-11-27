#include "Controller.h"

long int Controller::iOuter = 0;

Controller::Controller()
{
	Read_InputFile();
	if(SolverCode == 0)
	{
		eq[0] = 1, eq[1] = 1, eq[2] = 1;
		if(buoyancy == 1)
			M = new MomDriver(rho, mu, gx, gy, Beta);
		else
			M = new MomDriver(rho, mu);
		P = new PresDriver(rho);
		E = new ScalDriver(rho,mu,Pr,Cp); 
	}

	if(SolverCode == 1)
	{
		eq[0] = 1, eq[1] = 1;
		M = new MomDriver(rho, mu);
		P = new PresDriver(rho);
	}

	if(SolverCode == 2)
	{
		eq[0] = 1;
		E = new ScalDriver(rho,mu,Pr,Cp);
	}

	G = new Geometry();
}

Controller::~Controller()
{
	if(M) delete M;
	if(P) delete P;
	if(E) delete E;

	rho = 0.0, mu = 0.0, gx = 0.0, Beta = 0.0, Cp = 0.0, Pr = 0.0;
	Time_Scheme = 0.0;

	f.close();
	delete G;
	G = NULL;
}

void Controller::Read_InputFile()
{
	fstream g;
	string s;
	g.open("Flow_input.dat",ios::in);
	getline(g,s);
	g>>SolverCode;
	g.seekg(1, ios::cur);

	if(SolverCode == 0)	Neq = 3;
	else if(SolverCode == 1) Neq = 2;
	else if(SolverCode == 2) Neq = 1;

	eq = new unsigned short[Neq];
	Tolerance = new double[Neq];
	Inners = new long int[Neq];
	Lambda = new double[Neq];
	Urf = new double[Neq];

	getline(g,s);
	g>>rho>>mu>>Cp>>gx>>gy>>Beta>>Pr>>buoyancy;	g.seekg(1, ios::cur);
	getline(g,s);
	g>>Time_Scheme>>RestartFlag>>Ntime>>DataIntervals>>dt; g.seekg(1, ios::cur);
	getline(g,s);
	g>>MaxIterations>>StoreResiduals>>Coupling; g.seekg(1, ios::cur);
	getline(g,s);
	for(unsigned short t = 0; t < Neq; t++)
		g>>Inners[t]>>Tolerance[t]>>Lambda[t]>>Urf[t];

	g.close();
}

void Controller::Initialize()
{
	Driver::SetTimeScheme(Time_Scheme);
	Driver::Setrestart(RestartFlag);
	Driver::SetNtime(Ntime);
	Driver::SetTStep(dt);

	G->Initialize();

	if(SolverCode == 1 || SolverCode == 0)
	{
		M->Initialize(G);
		P->Initialize(G);

		P->SetPGrad(M->GetPGrad(X),X);P->SetPGrad(M->GetPGrad(Y),Y);
		P->SetVel(M->GetPhi(X),X);P->SetVel(M->GetPhi(Y),Y);
		P->SetVelC(M->GetVel(X),X);P->SetVelC(M->GetVel(Y),Y);
		P->SetAp(M->GetDiag(X),X);P->SetAp(M->GetDiag(Y),Y);

		M->SetPressure(P->GetPressure());
		M->SetmdotE(P->GetmdotE());
		M->SetmdotN(P->GetmdotN());

		P->Init();

		if(SolverCode == 0)
		{
			E->Initialize(G);
			if(buoyancy == 1)	M->SetFVar(E->GetPhi(0));

			E->SetmdotE(P->GetmdotE());
			E->SetmdotN(P->GetmdotN());	
			E->Init();
		}		
	}

	if(StoreResiduals == 1)
		f.open("Residuals.dat",ios::out);
}

bool Controller::PrintResiduals()
{
	long int I = Controller::GetOuter();
	if(I == 0)
	{
		cout<<"			Residuals"<<endl;
		if(SolverCode == 0)
			cout<<"-------Iter--------U----------V----------P----------T---------"<<endl;
		else if(SolverCode == 1)
			cout<<"-------Iter--------U----------V----------P----------"<<endl;
		return false;
	}
	if(I%100 == 0)
	{
		if(SolverCode == 0)
			cout<<"-------Iter--------U----------V----------P----------T---------"<<endl;
		else if(SolverCode == 1)
			cout<<"-------Iter--------U----------V----------P----------"<<endl;
	}

	double R = 0.0;
	if(SolverCode == 0 || SolverCode == 1)
	{
		double RT = 0.0;
		double Ru = M->GetRes(X); 
		double Rv = M->GetRes(Y);
		double Rp = P->GetRes();
		cout<<I<<"	"<<Ru<<"	"<<Rv<<"	"<<Rp;
		if(StoreResiduals == 1)	f<<I<<"	"<<Ru<<"	"<<Rv<<"	"<<Rp;
		R = max(Ru,Rv); R = max(R,Rp); 
		if(SolverCode == 0)	
		{
			double RT = E->GetRes();
			cout<<"	"<<RT;
			if(StoreResiduals == 1)	f<<"	"<<RT;
			R = max(R,RT);
		}
		cout<<endl;

		if(StoreResiduals == 1) f<<endl;
	}

	if(R < Coupling)
		return true;
	return false;
}

void Controller::Run()
{
	PrintResiduals();
	if(Time_Scheme == 2)
	{
		while(Controller::GetOuter() < MaxIterations)
		{
			if(SolverCode == 0)	E->Solve();
			M->Solve();
			P->Solve();
			Controller::IncOuter();
			if(PrintResiduals())
				break;
		}
	}
	else
	{
		long int itime = 0;
		double time;
		while(itime < Ntime)
		{
			while(Controller::GetOuter() < MaxIterations)
			{
				if(SolverCode == 0)	E->Solve();
				M->Solve();
				P->Solve();
				Controller::IncOuter();
				if(PrintResiduals())
					break;
			}

			if(itime%DataIntervals == 0)
			{
				M->Store();
				P->Store();
				E->Store();
			}

			itime++;
			time = time+dt;
		}
	}
}

int main()
{
	// Controller *C = new Controller();
	// C->Initialize();
	// C->Run();
	// if(C) delete C;

	Geometry *G = new Geometry();
	G->Initialize();

	// MomDriver *M = new MomDriver(1.0, 0.01, 0.0, 1.0, 100.00);
	// M->Initialize(G);

	MomDriver *M = new MomDriver(1.0, 0.001);
	M->Initialize(G);

	PresDriver *P = new PresDriver(1.0);
	P->Initialize(G);

	// ScalDriver *E = new ScalDriver(1.0, 0.01, 0.1);
	// E->Initialize(G);

	P->SetPGrad(M->GetPGrad(X),X);P->SetPGrad(M->GetPGrad(Y),Y);
	P->SetVel(M->GetPhi(X),X);P->SetVel(M->GetPhi(Y),Y);
	P->SetVelC(M->GetVel(X),X);P->SetVelC(M->GetVel(Y),Y);
	P->SetAp(M->GetDiag(X),X);P->SetAp(M->GetDiag(Y),Y);

	M->SetPressure(P->GetPressure());
	M->SetmdotE(P->GetmdotE());
	M->SetmdotN(P->GetmdotN());
	// M->SetFVar(E->GetPhi(0));

	// E->SetmdotE(P->GetmdotE());
	// E->SetmdotN(P->GetmdotN());

	Driver::Setrestart(1);
	P->Init();
	// E->Init();

	fstream f;
	f.open("Max_Residuals.dat",ios::out);
	double Ru = 0.0;
	int I = 1;
	cout<<"			Residuals"<<endl;
	cout<<"-------Iter--------U----------V----------P----------T---------"<<endl;
	while(1)
	{
		// E->Solve();
		M->Solve();
		P->Solve();
		Solver::im();
		double Ru = M->GetRes(X); 
		double Rv = M->GetRes(Y);
		double Rp = P->GetRes();
		// double RT = E->GetRes();
		cout<<I<<"	"<<Ru<<"	"<<Rv<<"	"<<Rp;
		// cout<<"	"<<RT;
		cout<<endl;
		
		f<<I<<"	"<<Ru<<"	"<<Rv<<"	"<<Rp<<endl;
		I++;
			
		// cout<<endl;
		Ru = max(Ru,Rv); Ru = max(Ru,Rp); 
		// Ru = max(Ru,RT);
		if(I %100 == 0)
			cout<<"-------Iter--------U----------V----------P----------T---------"<<endl;
		if(Ru < 1e-5)
		{
			cout<<"converged"<<endl;
			break;
		}

	}
	f.close();

	M->Store();
	P->Store();
	// E->Store();

	// if(E) delete E;
	// E = NULL;
	if(P) delete P;
	P = NULL;
	if(M) delete M;
	M = NULL;
	if(G) delete G;
	G = NULL;
}
