/*defines the simulation domain*/
#include <random>
#include <math.h>
#include "World.h"
#include "Field.h"
#include "Species.h"
	
//make an instance of the Rnd class
Rnd rnd;

using namespace std;

/*constructor*/
World::World(int ni, int nj, int nk):
	ni{ni}, nj{nj}, nk{nk},	nn{ni,nj,nk},
	phi(ni,nj,nk),rho(ni,nj,nk),node_vol(ni,nk,nk),
	ef(ni,nj,nk), object_id(ni,nj,nk)	{
		time_start =  chrono::high_resolution_clock::now();	//save starting time point
	}

/*sets domain bounding box and computes mesh spacing*/
void World::setExtents(double3 _x0, double3 _xm) {
	/*set origin and the opposite corner*/
	x0 = _x0;
	xm = _xm;

	/*compute spacing by dividing length by the number of cells*/
	for (int i=0;i<3;i++)
		dh[i] = (xm(i)-x0(i))/(nn(i)-1);

	//compute centroid
	xc = 0.5*(x0+xm);

	/*recompute node volumes*/
	computeNodeVolumes();
}

/*returns elapsed wall time in seconds*/
double World::getWallTime() {
  auto time_now = chrono::high_resolution_clock::now();
  chrono::duration<double> time_delta = time_now-time_start;
  return time_delta.count();
}

/*computes charge density from rho = sum(charge*den)*/
void World::computeChargeDensity(vector<Species> &species)
{
	rho = 0;
	for (Species &sp:species)
	{
		if (sp.charge==0) continue;	//don't bother with neutrals
		rho += sp.charge*sp.den;
	}
}	

/*computes node volumes, dx*dy*dz on internal nodes and fractional
 * values on domain boundary faces*/
void World::computeNodeVolumes() {
	for (int i=0;i<ni;i++)
		for (int j=0;j<nj;j++)
			for (int k=0;k<nk;k++)
			{
				double V = dh[0]*dh[1]*dh[2];	//default volume
				if (i==0 || i==ni-1) V*=0.5;	//reduce by two for each boundary index
				if (j==0 || j==nj-1) V*=0.5;
				if (k==0 || k==nk-1) V*=0.5;
				node_vol[i][j][k] = V;
			}
}

/* computes total potential energy from 0.5*eps0*sum(E^2)*/
double World::getPE() {
	double pe = 0;
	for (int i=0;i<ni;i++)
		for (int j=0;j<nj;j++)
			for (int k=0;k<nk;k++)
			{
				double3 efn = ef[i][j][k];	//ef at this node
				double ef2 = efn[0]*efn[0]+efn[1]*efn[1]+efn[2]*efn[2];

				pe += ef2*node_vol[i][j][k];
			}
	return 0.5*Const::EPS_0*pe;
}

/*sugarcubes a sphere centered at (x0,y0,z0)*/
void World::addSphere(double3 x0, double radius, double phi_sphere)
{
    /*save sphere parameters*/
    sphere_x0 = x0;
    sphere_rad2 = radius*radius;

    for (int i=0;i<ni;i++)
        for (int j=0;j<nj;j++)
            for (int k=0;k<nk;k++)
            {
                /*compute node position*/
                double3 x = pos(i,j,k);
                if (inSphere(x))
                {
                    object_id[i][j][k] = 1;
                    phi[i][j][k] = phi_sphere;
                }
            }
}

void World::addTubeBox(double3 min, double3 max, double holeRadius, double phi_object){
	
	for (int i=0;i<ni;i++)
        for (int j=0;j<nj;j++)
            for (int k=0;k<nk;k++)
            {
                /*compute node position*/
                double3 x = pos(i,j,k);
                if (inBox(x, min, max, holeRadius))
                {
                    object_id[i][j][k] = 1;
                    phi[i][j][k] = phi_object;
                }
            }
}

/*marks k=0 plane as 0V Dirichlet boundary*/
void World::addInlet() {
	for (int i=0;i<ni;i++)
		for (int j=0;j<nj;j++)
		{
			object_id[i][j][0] = 2;
			phi[i][j][0] = 0;
		}
}
	
/*returns true if point x is inside or on the sphere*/
bool World::inSphere(double3 x)
{
	double3 r = x-sphere_x0;	//ray to x

    double r_mag2 = (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    if (r_mag2<=sphere_rad2) return true;
    return false;
}

bool World::inBox(double3 point, double3 min, double3 max, double holeRadius) {
	//returns TRUE when in the box
	double r_mag2;
	double cyl_mag2;
	r_mag2 = (point[0] * point[0]) + (point[1] * point[1] );
	cyl_mag2 = (holeRadius * holeRadius);
	if ( point[0] >= min[0] && point[0] <= max[0] && 
		 point[1] >= min[1] && point[1] <= max[1] && 
		 point[2] >= min[2] && point[2] <= max[2] && 
		 r_mag2 >= cyl_mag2) return true;
	else return false;

}

bool World::inCylinder(double3 point, double holeRadius) {
	//this assumes the hole of radius R is centered at z = 0
	//TRUE when IN the cylinder
	double r_mag2;
	double cyl_mag2;
	r_mag2 = (point[0] * point[0]) + (point[1] * point[1] );
	cyl_mag2 = (holeRadius * holeRadius);
	if ( r_mag2 <= cyl_mag2) return true;
	else return false;
}
