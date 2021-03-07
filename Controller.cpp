#include "Controller.h"

long int Controller::iOuter = 0;

Controller::Controller()
{
	Read_InputFile();
	G = new Geometry();
	G->Initialize();
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
	Driver::Setrestart(RestartFlag);
	Driver::SetNtime(Ntime);
	Driver::SetTStep(dt);

	double iUrf[2], iLambda[2];

	if(SolverCode == 1 || SolverCode == 0)
	{
		iLambda[0] = Lambda[0], iLambda[1] = Lambda[0];
		iUrf[0] = Urf[0], iUrf[1] = Urf[0];
		M->SetTimeScheme(Time_Scheme);
		M->Initialize(G);
		M->SetParams(iLambda,iUrf,Tolerance[0],Inners[0]);

		iLambda[0] = Lambda[1];
		iUrf[0] = Urf[1];
		P->SetTimeScheme(Time_Scheme);
		P->Initialize(G);
		P->SetParams(iLambda,iUrf,Tolerance[1],Inners[1]);
		
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
			iLambda[0] = Lambda[2];
			iUrf[0] = Urf[2];
			E->SetTimeScheme(Time_Scheme);
			E->Initialize(G);
			E->SetParams(iLambda,iUrf,Tolerance[2],Inners[2]);

			if(buoyancy == 1)	
				M->SetFVar(E->GetPhi(0));

			E->SetmdotE(P->GetmdotE());
			E->SetmdotN(P->GetmdotN());	
			E->SetTimeScheme(Time_Scheme);
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
			if(SolverCode == 0)	
				E->Solve();
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
			if(SolverCode == 0)	
				E->Solve();
			M->Solve();
			P->Solve();
			Controller::IncOuter();
			if(PrintResiduals())
				break;
			
			if(itime%DataIntervals == 0)
			{
				M->Store();
				P->Store();
				if(SolverCode == 0)	E->Store();
			}

			itime++;
			time = time+dt;
			Driver::IncTime();
		}
	}

	if(SolverCode == 0) E->Store();
	M->Store();
	P->Store();
}

int main()
{
	Controller *C = new Controller();
	C->Initialize();
	C->Run();
	if(C) delete C;
}
