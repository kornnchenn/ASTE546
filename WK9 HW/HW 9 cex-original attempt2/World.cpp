/*defines the simulation domain*/
#include <random>
#include <math.h>
#include <iostream>
#include "World.h"
#include "Field.h"
#include "Species.h"
	
//make an instance of the Rnd class
Rnd rnd;

using namespace std;

/*constructor*/
World::World(int ni, int nj, int nk):
	ni{ni}, nj{nj}, nk{nk},	nn{ni,nj,nk},
	phi(nn),rho(nn),node_vol(nn),
	ef(nn), object_id(nn)	{
		time_start =  chrono::high_resolution_clock::now();	//save starting time point
	}

/*sets domain bounding box and computes mesh spacing*/
void World::setExtents(double3 _x0, double3 _xm) {
	/*set origin and the opposite corner*/
	x0 = _x0;
	xm = _xm;

	/*compute spacing by dividing length by the number of cells*/
	for (int i=0;i<3;i++)
		dh[i] = (xm[i]-x0[i])/(nn[i]-1);

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


void World::addBox(const double3 &box_x0, const double3 &box_xm, double phi_sphere){
	  	
	_box_x0 = box_x0;
	_box_xm = box_xm;

	for (int i=0;i<ni;i++)
        for (int j=0;j<nj;j++)
            for (int k=0;k<nk;k++)
            {
                /*compute node position*/
                double3 x = pos(i,j,k);
                if (inBox(x))
                {
                    object_id[i][j][k] = 1;
                    phi[i][j][k] = phi_sphere;
                }
            }

}

/*marks k=0 plane as 0V Dirichlet boundary*/
void World::addInlet(const double3 center, const double radius) {
	for (int i=0;i<ni;i++)
		for (int j=0;j<nj;j++)
		{
			//object_id[i][j][0] = 2;
			//phi[i][j][0] = 0;
			double dx = i - center[0];
			double dy = j - center[1];
			double distance = sqrt(dx * dx + dy * dy);

			if (distance <= radius) {
				object_id[i][j][0] = 2;
				phi[i][j][0] = 0;
			}

		}
}

bool World::inBox(const double3 &x) {
	//returns TRUE when in the box
    if ( x[0] >= _box_x0[0] && x[0] <= _box_xm[0] && 
         x[1] >= _box_x0[1] && x[1] <= _box_xm[1] && 
         x[2] >= _box_x0[2] && x[2] <= _box_xm[2]) return true;
    else return false;

}

/*checks for steady state by comparing change in mass, momentum, and energy*/
bool World::steadyState(vector<Species> &species) {
	// do not do anything if already at steady state
	if (steady_state) return true;

	double tot_mass = 0;
	double tot_mom = 0;
	double tot_en = getPE();
	for (Species &sp:species)
	{
		tot_mass += sp.getRealCount();	//number of real molecules
		double3 mom = sp.getMomentum();
		tot_mom += abs(mom[2]);		//z-component of momentum
		tot_en += sp.getKE();		//add kinetic energy
	}

	/*compute new values to last*/
	const double tol = 5e-4;
	if (abs((tot_mass-last_mass)/tot_mass)<tol &&
		abs((tot_mom-last_mom)/tot_mom)<tol &&
		abs((tot_en-last_en)/tot_en)<tol) {
		steady_state = true;
		cout<<"Steady state reached at time step "<<ts<<endl;
	}

	/*update prior values*/
	last_mass = tot_mass;
	last_mom = tot_mom;
	last_en = tot_en;
	return steady_state;
}

