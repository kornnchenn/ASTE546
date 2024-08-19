/*
ASTE-546 Simple Vlasov solver example
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <cstring>
#include <cmath>
#include <vector>

using namespace std;

//storage of domain parameters
struct World
{
	double L;	//length of x
	double v_max;	//extends in the velocity space
	double dx,dv;	//cell spacing
	int ni,nj;	//number of nodes in x,y
	
	void setLimits(double L, double v_max) {this->L=L;this->v_max=v_max;}
	void setNodes(int N, int M) {ni=N; nj=2*M-1;dx=L/(ni-1);dv=2*v_max/(nj-1);} 
	double getX(int i) {return 0+i*dx;}
	double getV(int j) {return -v_max+j*dv;}
	
	//linear interpolation: a higher order scheme needed!
	double interp(double **f, double x, double v)
	{
		//this version returns zero if out of bounds
		double fi = (x-0)/dx;
		double fj = (v-(-v_max))/dv;		
		if (fi<0 || fi>ni-1) return 0;
		if (fj<0 || fj>nj-1) return 0;
		
		int i = (int)fi;
		int j = (int)fj;
		double di = fi-i;
		double dj = fj-j;
		double val = (1-di)*(1-dj)*f[i][j];
		if (i<ni-1) val+=(di)*(1-dj)*f[i+1][j];
		if (j<nj-1) val+=(1-di)*(dj)*f[i][j+1];
		if (i<ni-1 && j<nj-1) val+=(di)*(dj)*f[i+1][j+1];
		return val;
		
	}
	
	~World() {}
};

//filter to eliminate garbage values like 1e-120
double filter(double a) {if (std::abs(a)<1e-20) return 0; else return a;}

/*saves the provided scalars and vectors to a file*/
void saveVTK(int time_step, World &world, map<string,double**> scalars2D, map<string,double*> scalars1D)
{
	//generate file name
	stringstream ss;
	ss<<"results/vlasov";
	if (time_step>=0)
		ss<<"_"<<setw(6)<<setfill('0')<<time_step;
	ss<<".vti";	
	ofstream out(ss.str());

	out<<setprecision(4);

	out<<"<VTKFile type=\"ImageData\">\n";
	out<<"<ImageData Origin=\""<<0<<" "<<-world.v_max<<" "<<0<<"\"";
	out<<" Spacing=\""<<world.dx<<" "<<world.dv<<" "<<1<<"\"";
	out<<" WholeExtent=\""<<0<<" "<<world.ni-1<<" "<<0<<" "<<world.nj-1<<" "<<0<<" "<<0<<"\">\n";
	out<<"<Piece Extent=\""<<0<<" "<<world.ni-1<<" "<<0<<" "<<world.nj-1<<" "<<0<<" "<<0<<"\">\n";
	out<<"<PointData>\n";
		
	//user vars, p.first is the string key, p.second is the double* pointer to data
	for (pair<string,double**> p : scalars2D)
	{
		//p.first is the string key, p.second is the double* pointer to data
		out<<"<DataArray Name=\""<<p.first<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		
		for (int j=0;j<world.nj;j++)
		{
			for (int i=0;i<world.ni;i++)
			{
				out<<setprecision(4)<<filter(p.second[i][j])<<" ";					
			}
			out<<"\n";
		}
		out<<"</DataArray>\n";
	}

	for (pair<string,double*> p : scalars1D)
	{
		//p.first is the string key, p.second is the double* pointer to data
		out<<"<DataArray Name=\""<<p.first<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		
		for (int j=0;j<world.nj;j++)
		{
			for (int i=0;i<world.ni;i++)
			{
				out<<setprecision(4)<<filter(p.second[i])<<" ";					
			}
			out<<"\n";
		}
		out<<"</DataArray>\n";
	}

	out<<"</PointData>\n";
	out<<"</Piece>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";

	out.close();
}

/*global container to store all allocated memory*/
vector<pair<int,double**>> _vecsAllocated;

/*allocates a 2D [ni][nj] array*/
double** newAndClear(int ni, int nj)
{
	//first allocate a 1D array for the first index
	double **p = new double*[ni];
	for (int i=0;i<ni;i++)
	{
		p[i] = new double[nj];
		memset(p[i],0,sizeof(double)*nj);	
	}
	
	//add to container
	_vecsAllocated.emplace_back(pair<int,double**>(ni,p));
	
	return p;
}

/*allocates a 1D [ni] array*/
double* newAndClear(int ni)
{
	//first allocate a 1D array for the first index
	double **p = new double*[1];
	p[0] = new double[ni];
	memset(p[0],0,sizeof(double)*ni);	
	
	//add to container
	_vecsAllocated.emplace_back(pair<int,double**>(1,p));
	
	return p[0];
}


/*deletes all memory allocated by newAndClear, gets pointers and sizes from _vecsAllocated*/
void deleteAll()
{
	//for-each loop from C++11, identical to: for (int i=0;i<v.size();i++) {auto &a=v[i];...}
	for (auto &p : _vecsAllocated)
	{
		//delete the array of pointers
		for (int i=0;i<p.first;i++)
			delete[] p.second[i];
		delete[] p.second;
	}
}

/*solves Poisson's equation with Dirichlet boundaries using the direct Thomas algorithm
and returns the electric field*/
void solvePoissonsEquation(World &world, double *ne, double *E)
{
	double *a = new double[world.ni];
	double *b = new double[world.ni];
	double *c = new double[world.ni];
	double *d = new double[world.ni];
	double *x = new double[world.ni];
	int ni = world.ni;
	
	for (int i=0;i<ni;i++)
	{
		if (i==0 || i==ni-1)
		{
			//dirichlet boundary
			b[i] = 1;
			d[i] = (i==0)?15:0;		//0V on left, 5 on right
		}
		else
		{
			double dx2 = world.dx*world.dx;
			a[i] = 1/dx2;
			b[i] = -2/dx2;
			c[i] = 1/dx2;
			d[i] = 1-ne[i];		//rho
		}
	}
	
	//initialize
	c[0] = c[0]/b[0];
	d[0] = d[0]/b[0];
	
	//forward
	for (int i=1;i<ni;i++)
	{
		if (i<ni-1)
			c[i] = c[i]/(b[i]-a[i]*c[i-1]);
		
		d[i] = (d[i] - a[i]*d[i-1])/(b[i] - a[i]*c[i-1]);
	}
	
	//backward
	x[ni-1] = d[ni-1];
	for (int i=ni-2;i>=0;i--)
	{
		x[i] = d[i] - c[i]*x[i+1];
	}
		
	//set E field
	E[0] = -(x[1]-x[0])/world.dx;
	E[ni-1] = -(x[ni-1]-x[ni-2])/world.dx;
	for (int i=1;i<ni-1;i++)
	{
		E[i] = -(x[i+1]-x[i-1])/(2*world.dx);
	}
	
	delete[] a;
	delete[] b;
	delete[] c;
	delete[] d;
	delete[] x;
}


int main()
{
	//constants and parameters
	const double pi = acos(-1.0);		//pi

	//create a variable of type World
	World world;	
	world.setLimits(10,10);
	world.setNodes(21,21);
	int ni = world.ni;
	int nj = world.nj;
	double dx = world.dx;
	double dv = world.dv;
	double dt = 1/8.0;

	
	cout<<"dx: "<<dx<<" "<<"dv: "<<dv<<endl;
			
	double **f = newAndClear(ni,nj); //f
	double **fs = newAndClear(ni,nj); //f
	double **fss = newAndClear(ni,nj); //f
	double *ne = newAndClear(ni);	//number density, 1D vector
	double *E = newAndClear(ni);	//number density, 1D vector
	
	//map is a list of keys and corresponding values
	map<string,double**> scalars2D; 
	map<string,double*> scalars1D; 
	
	scalars2D["f"] = f;	
	scalars1D["ne"] = ne;
	scalars1D["E"] = E;
	
	//set initial distribution
	for (int i=0;i<ni;i++)
		for (int j=0;j<nj;j++)
		{
			double x = world.getX(i);
			double v = world.getV(j);
		
			double f0 = 1.0/sqrt(2*pi)*exp(-v*v/2.0);
			f[i][j] = f0*(exp(-(i-0.5*ni)*(i-0.5*ni)/2.0));		
			f[i][j] = f0;
			
		}
	
	//set some constant e field
	for (int i=0;i<ni;i++)
		E[i] = 0;
	
	int it;
	
	for (it=0;it<100;it++)
	{
		if (it%2==0)
			saveVTK(it,world,scalars2D,scalars1D);

		//compute f*
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
			{
				double v = world.getV(j);
				double x = world.getX(i);
				
				fs[i][j] = world.interp(f,x-v*0.5*dt,v);								
			}
			
		//compute number density by integrating f with the trapezoidal rule		
		for (int i=0;i<ni;i++)
		{
			ne[i] = 0;
			for (int j=0;j<nj-1;j++)
				ne[i]+=0.5*(fs[i][j+1]+fs[i][j])*dv;
		}

		
		//different schemes to obtain the E field
		
		//forward difference
		/*
		E[0] = 0;
		for (int i=1;i<ni;i++)
		{
			E[i] = E[i-1] + dx*(1-ne[i-1]);
		}
		*/
		
		//backward difference
		/*
		E[ni-1] = 0;
		for (int i=ni-2;i>=0;i--)
		{
			E[i] = E[i+1] + dx*(1-ne[i+1]);
		}
		*/
		
		//solution of the Poisson's equation
		solvePoissonsEquation(world,ne,E);
		
		//compute f**
		for (int i=0;i<ni;i++)
			for(int j=0;j<nj;j++)
			{
				double v = world.getV(j);
				double x = world.getX(i);
				fss[i][j] = world.interp(fs,x,v+E[i]*dt);				
			}
		
		//compute f(n+1)
		for (int i=0;i<ni;i++)
			for(int j=0;j<nj;j++)
			{
				double v = world.getV(j);
				double x = world.getX(i);
				f[i][j] = world.interp(fss,x-v*0.5*dt,v);
			}		
			
	}
		
		
	saveVTK(it,world,scalars2D,scalars1D);
		
	return 0;
}
