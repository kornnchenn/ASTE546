//sheath 1d
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <sstream>
#include <iomanip>

using namespace std;
using dvector = vector<double>;

// specify new data object of type Rnd
class Rnd {
public:
    // overload the () operator to return random value
	double operator() () {
		return dist(gen);
	}

protected:
 // class member variables, not accessible from outside the class
 std::default_random_engine gen{0*std::random_device()()};  // must use {} here
 std::uniform_real_distribution<double> dist{0.0,1.0};
};

Rnd rnd;      // create a variable of type Rnd

struct World {

	// convert physical position to grid logical coordinate
	double XtoL(double x) {
		return (x-x0)/dx;
	}

	double x0;    // mesh origin
	double xm;    // maximum point
	double dx;    // cell size;
	int ni;       // number of mesh nodes

	int ts;       // current time step
	double dt;    // time step size (seconds)
};

struct Particle {
	double x, v;
};

class Species {
public:
	// constructor
	//Species() {}

	Species(double m, double q, double mpw, size_t N) {
		(*this).m = m;
		 this->q = q;
		 this->mpw = mpw;
		 np = 0;
		 np_alloc = N;
		 part = new Particle[np_alloc];
	}


	// destructor
	~Species() {
		delete[] part;
		part = nullptr;
	}

  double m, q;
  double mpw;

  Particle* operator[](int i) {
	  if (i>=0 && i<np_alloc)
		  return &part[i];
	  else
		  throw runtime_error("out of bounds!");
  }

//protected:
  size_t np_alloc;
  size_t np;
  Particle *part;

};



/*function prototypes*/
//bool outputCSV(const dvector &phi, const dvector &rho, const dvector &ef,
//	  const dvector &ndi, const dvector &nde, World &world);
void solvePotentialDirect(dvector &phi, const dvector &rho, World &world);
bool solvePotentialGS(dvector &phi, const dvector &rho, World &world, int max_it);
void computeEF(dvector &ef, const dvector &phi, World &world, bool second_order);
void advanceSpecies(Species &species, dvector &ef, World &world);
void rewindVelocity(Species &species,  dvector &ef, World &world);
void computeNumberDensity(Species &species, dvector &nd, World &world);


double gather(double li, const dvector &f) {
	int i = (int) li;
	double di = li-i;
	return f[i]*(1-di) + f[i+1]*di;
}

void scatter(double li, dvector &f, double val) {
	int i = (int) li;
	double di = li-i;
	f[i]+=(1-di)*val;
	f[i+1]+=di*val;
}

/*constants*/
namespace Const
{
	const double QE = 1.602176565e-19;	// C, electron charge
	const double EPS_0 = 8.85418782e-12;// C/V/m, vacuum permittivity
	const double ME = 9.10938215e-31;	// kg, electron mass
	const double AMU = 1.660539e-27;    // mass of a proton
};

using namespace Const;	//to avoid having to write Const::QE



/*main*/
int main()
{
	World world;

	world.ni = 41;	//number of nodes
	world.x0 = 0;	//origin
	world.xm = 0.1;	//opposite end
	world.dx = (world.xm-world.x0)/(world.ni-1);	//node spacing

	vector<double> phi(world.ni);
	dvector rho(world.ni,QE*1e12);
	dvector ef(world.ni);

	constexpr int N = 100000; //number of particles?

	constexpr double n0 = 1e15;    // the desired number density

	double n_real = n0*(world.xm-world.x0);   // number of real particles in domain of volume (xd-x0)*1*1

	Species ions = {16*Const::AMU, Const::QE, n_real/N, N};
	Species eles(Const::ME, -Const::QE, n_real/N, N);

   double vth_i = 500;
	
    // inject stationary particles
    for (int p=0;p<ions.np_alloc;p++) {

    	Particle *part = ions[p];

    	part->x = world.x0 + rnd()*(world.xm-world.x0);
    	part->v = 0;  //stationary

    	part->v = vth_i*(rnd()+rnd()+rnd()-1.5);
    	
    	ions.np++;   //increment counter of particles
    }

	double vth_e = 5e5;
    // inject stationary particles
    for (int p=0;p<eles.np_alloc;p++) {
    	eles.part[p].x = world.x0 + rnd()*(world.xm-world.x0);
    	eles.part[p].v = 0;  //stationary
    	eles.part[p].v = vth_e*(rnd()+rnd()+rnd()-1.5);  
    	
    	eles.np++;
    }

	world.dt = 1e-10;//2e-8; //5e-11;


	vector<double> ndi(world.ni);
	dvector nde(world.ni);

	/*rewind for leapfrog!*/

	advanceSpecies(ions,ef,world);
	advanceSpecies(eles,ef,world);
	computeNumberDensity(ions,ndi,world);
	computeNumberDensity(eles,nde,world);

	for (int i=0;i<world.ni;i++) {
		rho[i] = ndi[i]*ions.q + nde[i]*eles.q;
	}
	
	solvePotentialDirect(phi, rho, world);
	computeEF(ef, phi, world, true);
	rewindVelocity(ions,ef,world);
	rewindVelocity(eles,ef,world);
	
	// simulation main loop
    ofstream out("diag.csv");
    out << "ts, time, num_ions, num_electrons" << endl;
    double sim_time;

	for (world.ts=0;world.ts<20000;world.ts++) {

		advanceSpecies(ions,ef,world);
		advanceSpecies(eles,ef,world);

		computeNumberDensity(ions,ndi,world);
		computeNumberDensity(eles,nde,world);

		for (int i=0;i<world.ni;i++) {
			rho[i] = ndi[i]*ions.q + nde[i]*eles.q;
		}

		solvePotentialDirect(phi, rho, world);
		computeEF(ef, phi, world, true);


		if (world.ts%5==0) {
            sim_time = world.ts * world.dt;
			out << world.ts <<","<< sim_time << "," << ions.np <<"," << eles.np << endl;
		}
		
		//if (world.ts%100==0)
			//ouput to a CSV file for plotting
			//outputCSV(phi,rho,ef,ndi,nde,world);
	}

	//cout<<"Done!"<<endl;


	return 0;	//normal exit
}

void advanceSpecies(Species &species, dvector &ef, World &world) {

	for (int p=0;p<species.np;p++) {
		double li = world.XtoL(species.part[p].x);       // get particle's grid index

		double E_p = gather(li, ef);            // electric field at particle position

		// here we go from -0.5 to +0.5
		species.part[p].v += species.q*E_p/species.m*world.dt;          // integrate velocity through dt


		// now v is at 0.5, use v at 0.5 to integrate x from 0 to 1
		species.part[p].x += species.part[p].v*world.dt;               // integrate position
	}
	
	/*perform particle removal sweep*/
	for (int p=0;p<species.np;p++) {
		
		if (species.part[p].x<world.x0 || species.part[p].x>=world.xm) {             // domain extent is [x0, xd)
			// copy particle from the end of the array to position p
			species.part[p] = species.part[species.np-1];    // this copies all struct members
			species.np--;		// reduce array size
			p--;				// reduce array counter to recheck part[p]
		}
	}
}

void rewindVelocity(Species &species, dvector &ef, World &world) {
	for (int p=0;p<species.np;p++) {
		double li = world.XtoL(species.part[p].x);       // get particle's grid index
		double E_p = gather(li, ef);            // electric field at particle position

		// here we go from -0.5 to +0.5
		species.part[p].v -= 0.5*species.q*E_p/species.m*world.dt;          // integrate velocity through dt
	}
}

void computeNumberDensity(Species &species, dvector &nd, World &world) {
	for (int i=0;i<nd.size();i++) nd[i] = 0;

	for (int p=0;p<species.np;p++) {
		double li = world.XtoL(species.part[p].x);
		scatter(li,nd,species.mpw);
	}
	
	size_t ni = nd.size();
	for (int i=0;i<ni;i++)
		nd[i] /= world.dx;
		
	nd[0] *= 2.0;
	nd[ni-1] *= 2.0;	
}

/*outputs the given fields to a CSV file, returns true if ok
bool outputCSV(const dvector &phi, const dvector &rho, const dvector &ef,
			   const dvector &ndi,const dvector &nde, World &world)
{
	stringstream name;
	name<<"results/field_"<<setw(6)<<setfill('0')<<world.ts<<".csv";
	
	ofstream out(name.str());	//open file for writing
	
	if (!out)
	{
		cerr<<"Could not open output file!"<<endl;
		return false;
	}

	out<<"x,phi,rho,ef,ndi,nde\n";		//write header
	for (int i=0;i<phi.size();i++)
	{
		out<<world.x0+i*world.dx; //write i-th position
		out<<","<<phi[i]<<","<< rho[i]<<","<<ef[i]<<","<<ndi[i]<<","<<nde[i]; //write values
		out<<"\n";	//new line, not using endl to avoid buffer flush
	}

	//file closed automatically when "out" variable is destroyed
	return true;
}
*/

/*solves Poisson's equation with Dirichlet boundaries using the Thomas algorithm*/

void solvePotentialDirect(dvector &phi, const dvector &rho, World &world)
{
	int ni = phi.size();	//number of mesh nodes
	dvector a(ni);			//allocate memory for the matrix coefficients
	dvector b(ni);
	dvector c(ni);
	dvector d(ni);

	//set coefficients
	for (int i=0;i<ni;i++)
	{
		if (i==0 || i==ni-1)	//Dirichlet boundary
		{
			b[i] = 1;		//1 on the diagonal
			d[i] = 0;		//0 V
		}
		else
		{
			a[i] = 1/(world.dx*world.dx);
			b[i] = -2/(world.dx*world.dx);
			c[i] = 1/(world.dx*world.dx);
			d[i] = -rho[i]/EPS_0;
		}
	}

	//initialize
	c[0] = c[0]/b[0];
	d[0] = d[0]/b[0];

	//forward step
	for (int i=1;i<ni;i++)
	{
		if (i<ni-1)
			c[i] = c[i]/(b[i]-a[i]*c[i-1]);

		d[i] = (d[i] - a[i]*d[i-1])/(b[i] - a[i]*c[i-1]);
	}

	//backward substitution
	phi[ni-1] = d[ni-1];
	for (int i=ni-2;i>=0;i--)
	{
		phi[i] = d[i] - c[i]*phi[i+1];
	}
}

/* solves potential using the Gauss Seidel Method, returns true if converged*/
bool solvePotentialGS(dvector &phi, const dvector &rho, World &world, int max_it)
{
	double L2;
	double dx2 = world.dx*world.dx;		//precompute dx*dx
	const double w = 1.4;
	int ni = phi.size();	//number of mesh nodes

	/*solve potential*/
	for (int solver_it=0;solver_it<max_it;solver_it++)
	{
		phi[0] = 0;			//dirichlet boundary on left
		phi[ni-1] = 0;		//dirichlet boundary on right

		/*Gauss Seidel method, phi[i-1]-2*phi[i]+phi[i+1] = -dx^2*rho[i]/eps_0*/
		for (int i=1;i<ni-1;i++)
		{
			double g = 0.5*(phi[i-1] + phi[i+1] + dx2*rho[i]/EPS_0);
			phi[i] = phi[i] + w*(g-phi[i]);	//SOR
		}

		/*check for convergence*/
		if (solver_it%50==0)
		{
			double sum = 0;

			//internal nodes, automatically satisfied on Dirichlet boundaries
			for (int i=1;i<ni-1;i++)
			{
				double R = -rho[i]/EPS_0 - (phi[i-1] - 2*phi[i] + phi[i+1])/dx2;
				sum+=R*R;
			}
			L2 = sqrt(sum/ni);
			if (L2<1e-6)
			{
				//cout<<"Gauss-Seidel converged after "<<solver_it<<" iterations"<<endl;
				return false;
			}
		}
	}
	//cout<<"Gauss-Seidel failed to converge, L2="<<L2<<endl;
	return true;
}
/* computes electric field by differentiating potential*/

void computeEF(dvector &ef, const dvector &phi, World &world, bool second_order)
{
	int ni = phi.size();	//number of mesh nodes

	//central difference on internal nodes
	for (int i=1;i<ni-1;i++)
		ef[i] = -(phi[i+1]-phi[i-1])/(2*world.dx);

	//boundaries
	if (second_order)
	{
		ef[0] = (3*phi[0]-4*phi[1]+phi[2])/(2*world.dx);
		ef[ni-1] = (-phi[ni-3]+4*phi[ni-2]-3*phi[ni-1])/(2*world.dx);
	}
	else	//first order
	{
		ef[0] = (phi[0]-phi[1])/world.dx;
		ef[ni-1] = (phi[ni-2]-phi[ni-1])/world.dx;
	}
}
