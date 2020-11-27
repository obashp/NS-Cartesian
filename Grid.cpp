#include "Grid.h"

Geometry::Geometry()
{
	x = NULL; y = NULL;
	X = NULL; Y = NULL;
	Nx = 0, Ny = 0;
	CNx = 0, CNy = 0;
	N = 0;
	Map21 = NULL;
	gx = NULL, gy = NULL;
}

void Geometry::Initialize()
{
	//Read the grid points from a file
	fstream f;
	string data;
	f.open("Grid.dat",ios::in);
	f>>Nx>>Ny;
	CNx = Nx+1; CNy = Ny+1;
	N = CNx*CNy;
	
	Map21 = new long int[CNx];
	x = new double[Nx]; y = new double[Ny];
	long int i = 0, j = 0;
	while(1)
	{
		f>>x[i++];
		if(i == Nx)
			break;
	}

	while(1)
	{
		f>>y[j++];
		if(j == Ny)
			break;
	}
	f.close();

	for(i = 0; i < CNx; i++)
		Map21[i] = i*CNy;

	X = new double[CNx]; Y = new double[CNy];
	gx = new double[CNx]; gy = new double[CNy];

	X[0] = x[0];
	gx[0] = 0.0;
	for(long int I = 1; I < Nx; I++)
	{
		X[I] = 0.5*(x[I] + x[I-1]);
		gx[I] = (I != Nx-1 ? (x[I]-x[I-1])/(x[I+1]-x[I-1]) : 1.0);
	}
	gx[CNx] = 0.0;
	X[CNx-1] = x[Nx-1];

	Y[0] = y[0];
	gy[0] = 0.0;
	for(long int J = 1; J < Ny; J++)
	{
		Y[J] = 0.5*(y[J] + y[J-1]);
		gy[J] = (J < Ny-1 ? (y[J]-y[J-1])/(y[J+1]-y[J-1]) : 1.0);
	}
	gy[CNy] = 0.0;
	Y[CNy-1] = y[Ny-1];

	

	f.open("Gcheck.dat",ios::out);
	for(long int I = 0; I < CNx; I++)
	{
		for(long int J = 0; J < CNy; J++)
			f<<Map21[I]+J<<"	"<<X[I]<<"	"<<Y[J]<<endl;
		f<<endl;
	}
	f.close();
}

Geometry::~Geometry()
{
	delete[] x,y, X,Y, gx, gy;
	x = NULL; y = NULL;
	X = NULL; Y = NULL;
	gx = NULL; gy = NULL;

	delete[] Map21;
	Map21 = NULL;

	Nx = 0, Ny = 0;
	CNx = 0, CNy = 0;
	N = 0;
	cout<<"Deleted Geometry object"<<endl;
}