#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <algorithm>

#ifndef GRID_H
#define GRID_H

using namespace std;

class Geometry
{
public:
	long int Nx, Ny;	// No. of nodes along each direction
	long int CNx, CNy;	// No. of cells along each direction. 
	long int N;	// No. of total cells
	long int *Map21;	//Map from 2-D indices to 1-D indices

	double *x, *y;	//Co-ordinates of x and y of nodes making up domain
	double *X, *Y;	//Coordinates of centroids of cells. 
	double *gx, *gy;

	Geometry();
	~Geometry();
	void Initialize();
};

#endif