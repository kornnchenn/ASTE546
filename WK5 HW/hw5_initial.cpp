//sheath 1d
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <sstream>
#include <iomanip>
#include <cstring>

using namespace std;
using dvector = vector<double>;

/*constants*/
namespace Const
{
	const double QE = 1.602176565e-19;	// C, electron charge
	const double EPS_0 = 8.85418782e-12;// C/V/m, vacuum permittivity
	const double ME = 9.10938215e-31;	// kg, electron mass
	const double AMU = 1.660539e-27;    // mass of a proton
	const double Kb = 1.380649e-23;		// Boltzmann constant
	const double EvToK = QE/Kb;			// 1 eV to Kelvin, ~11604K
};

using namespace Const;	//to avoid having to write Const::QE



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

	// constructor
	World (int ni, double x0, double xm) : phi(ni), rho(ni), ef(ni) {
		this->ni = ni;	  // alternatively can use ni{ni} in the initializer list
		this->x0 = x0;
		this->xm = xm;
		this->dx = (xm-x0)/(ni-1);
	}
	
	void setTime(double dt) {ts = 0; this->dt = dt;}
	
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
	
	dvector phi;
	dvector rho;
	dvector ef;
};

struct Particle {
	double x, v;
};

class Species {
public:
	// constructor
	//Species() {}

	Species(double m, double q, double mpw, size_t N, World &world) : 
		den(world.ni), vel(world.ni), world{world} {
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
		  throw runtime_error("out of bounds! ");
  }
  
  // compute average kinetic energy
  double getAveKE() {
  	if (np==0) return 0;
  	double sum = 0;
  	for (int i=0; i<np; i++) sum+=part[i].v*part[i].v;
  
  	return 0.5*m*sum/np;   
  }
  
  double getCurrent(World &world) {
  	double time = world.dt*(world.ts-ts_last);
  	double dnp = (np_last - np);  	
  	
  	double current; 
  	if (np_last>0) current = q*dnp*mpw/time; else current = 0; 	
  	
  	np_last = np;
  	ts_last = world.ts;
  	
  	return current;
  }

  double sampleVel(double T) {
    double v_th = sqrt(2*Const::Kb*T/m);
    constexpr int M=6;
    double R_sum = 0;
    for (int i=0;i<M;i++) R_sum +=rnd();
    return v_th/sqrt(2.0) * (R_sum-M/2.0)/sqrt(M/12.0);
}
  
  void advanceSpecies();
  void rewindVelocity();
  void computeNumberDensity();
  void computeVelocity();
  void outputParticles(string prefix, int count);

//protected:
  size_t np_alloc;
  size_t np;
  Particle *part;
  
  dvector den;   // number density
  dvector vel;

  
protected:
  int np_last = -1;
  double ts_last = -1;
   World &world;
};


class Results {
 public:

  // constructor,  
  Results(World &world, int nj) : ni{world.ni}, nj{nj}, world{world}  {
  	// allocate buffers and clear all bytes to 0
  	phi = new double[ni*nj]; memset(phi,0,ni*nj*sizeof(double));
  	rho = new double[ni*nj]; memset(rho,0,ni*nj*sizeof(double));
  	ef = new double[ni*nj]; memset(ef,0,ni*nj*sizeof(double));
  	ndi = new double[ni*nj]; memset(ndi,0,ni*nj*sizeof(double));
  	nde = new double[ni*nj]; memset(nde,0,ni*nj*sizeof(double));
  	vi = new double[ni*nj]; memset(vi,0,ni*nj*sizeof(double));
  	ve = new double[ni*nj]; memset(ve,0,ni*nj*sizeof(double));
  }
  
  // destructor, frees memory
  ~Results() {
  	delete[] phi; phi = nullptr;
  	delete[] rho; rho = nullptr;
  	delete[] ef; ef = nullptr;
	delete[] ndi; ndi = nullptr;
  	delete[] nde; nde = nullptr;
	delete[] vi; vi = nullptr;
  	delete[] ve; ve = nullptr;
  }

  void addData(Species &ions, Species &eles);
  bool outputVTI();

protected:  
  World &world; 		// reference to world
  
  // data buffers
  double *phi = nullptr, *rho = nullptr, *ef = nullptr;
  double *ndi = nullptr, *nde = nullptr;
  double *vi = nullptr, *ve = nullptr;
  
  // dimensions
  int ni = 0, nj = 0;       // dimensions
  int j = 0;
};


/*function prototypes*/
bool outputCSV(const dvector &ndi, const dvector &nde, World &world);
void solvePotentialDirect(World &world);
bool solvePotentialGS(World &world, int max_it);
void computeEF(World &world, bool second_order);

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


/*main*/
int main()
{
	World world(41,0,0.1); 		// number of nodes
	world.setTime(1e-11);
	int max_ts = 20000;
	
	constexpr int N = 400000;

	constexpr double n0 = 1e14;    // the desired number density

	double n_real = n0*(world.xm-world.x0);   // number of real particles in domain of volume (xd-x0)*1*1

//Original
	Species ions = {16*Const::AMU, Const::QE, n_real/N, N, world};
	Species eles(Const::ME, -Const::QE, n_real/N, N, world);

//Attempt to change charge
	//Species ions = {16*Const::AMU, 0.1 * EvToK, n_real/N, N, world};
	//Species eles(Const::ME, -1.0 * EvToK , n_real/N, N, world);

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

	int add_skip = 20;

	/*rewind for leapfrog!*/

	ions.advanceSpecies();
	eles.advanceSpecies();
	ions.computeNumberDensity();
	eles.computeNumberDensity();

	for (int i=0;i<world.ni;i++) {
		world.rho[i] = ions.den[i]*ions.q + eles.den[i]*eles.q;
	}
	
	solvePotentialDirect(world);
	computeEF(world, true);
	ions.rewindVelocity();
	eles.rewindVelocity();
	
	// diag file
	ofstream diag("diag.csv");
	diag<<"ts,time,num_ions,num_eles,KE_ions,KE_eles,I_ions, I_eles\n";
	
	Results results(world, max_ts/add_skip);
	
	// simulation main loop
	for (world.ts=0;world.ts<max_ts;world.ts++) {

		ions.advanceSpecies();
		eles.advanceSpecies();

		ions.computeNumberDensity();
		eles.computeNumberDensity();

		for (int i=0;i<world.ni;i++) {
			world.rho[i] = ions.den[i]*ions.q + eles.den[i]*eles.q;
		}

		solvePotentialDirect(world);
		computeEF(world, true);

		if (world.ts%50==0) {
			cout<<"ts: "<<world.ts<<", num_ions: "<<ions.np<<", num_eles: "<<eles.np<<endl;
		}
		
		if (world.ts%100==0) {
			diag<<world.ts<<","<<world.ts*world.dt<<","<<ions.np<<","<<eles.np;
			diag<<","<<ions.getAveKE()/Const::QE<<","<<eles.getAveKE()/Const::QE;
			diag<<","<<ions.getCurrent(world)<<","<<-eles.getCurrent(world)<<"\n";			
		}
		
		if (world.ts%add_skip==0)
			// add data to the buffer
			results.addData(ions,eles);
			
		if (world.ts%1000==0) {
			ions.outputParticles("ions",50000);
			eles.outputParticles("eles",50000);
		} 
		
	}

	results.outputVTI();
	
	cout<<"Done!"<<endl;

	return 0;	//normal exit
}

void Species::advanceSpecies() {

	for (int p=0;p<np;p++) {
		double li = world.XtoL(part[p].x);       // get particle's grid index

		double E_p = gather(li, world.ef);            // electric field at particle position

		// here we go from -0.5 to +0.5
		part[p].v += q*E_p/m*world.dt;          // integrate velocity through dt


		// now v is at 0.5, use v at 0.5 to integrate x from 0 to 1
		part[p].x += part[p].v*world.dt;               // integrate position
	}
	
	/*perform particle removal sweep*/
	for (int p=0;p<np;p++) {
		
		if (part[p].x<world.x0 || part[p].x>=world.xm) {             // domain extent is [x0, xd)
			// copy particle from the end of the array to position p
			part[p] = part[np-1];    // this copies all struct members
			np--;		// reduce array size
			p--;				// reduce array counter to recheck part[p]
		}
	}
}


void Species::rewindVelocity() {
	for (int p=0;p<np;p++) {
		double li = world.XtoL(part[p].x);    // get particle's grid index
		double E_p = gather(li, world.ef);            // electric field at particle position

		// here we go from -0.5 to +0.5
		part[p].v -= 0.5*q*E_p/m*world.dt;          // integrate velocity through dt
	}
}


void Species::computeNumberDensity() {
	size_t ni = den.size();
	for (int i=0;i<ni;i++) den[i] = 0;

	for (int p=0;p<np;p++) {
		double li = world.XtoL(part[p].x);
		scatter(li,den,mpw);
	}
	
	for (int i=0;i<ni;i++)
		den[i] /= world.dx;
		
	den[0] *= 2.0;
	den[ni-1] *= 2.0;	
}

void Species::computeVelocity() {
	size_t ni = vel.size();
	dvector count(ni);			// local variable to track sum(mpw)
	
	for (int i=0;i<ni;i++) vel[i] = 0;

	for (int p=0;p<np;p++) {
		double li = world.XtoL(part[p].x);
		scatter(li,vel,mpw*part[p].v);   // accumulate mpw*vel
		scatter(li,count,mpw);		 	// accumulate mpw
	}
	
	for (int i=0;i<ni;i++)               // divide sum(mpw*vel) by sum(mpw)
		if (count[i]>0) vel[i] /= count[i];  // avoid div by zero
		else vel[i] = 0;		
}

void Species::outputParticles(string prefix, int count) {
	stringstream name;
	name<<"results/"<<prefix<<"_"<<setw(6)<<setfill('0')<<world.ts<<".csv";
	
	ofstream out(name.str());	//open file for writing
	out<<"x,v\n";
	
	for (int p=0;p<np;p++) {
		// randomly select particles with probability count/np
		if (rnd()< count/(double)np) {  // note the cast to double to avoid integer division
			out<<part[p].x<<","<<part[p].v<<"\n";		
		}	
	}
}

/*outputs the given fields to a CSV file, returns true if ok*/
bool outputCSV(Species &ions, Species &eles, World &world)
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
	for (int i=0;i<world.ni;i++)
	{
		out<<world.x0+i*world.dx; //write i-th position
		out<<","<<world.phi[i]<<","<<world.rho[i]<<","<<world.ef[i]<<","<<ions.den[i]<<","<<eles.den[i]; //write values
		out<<"\n";	//new line, not using endl to avoid buffer flush
	}

	//file closed automatically when "out" variable is destroyed
	return true;
}

/*solves Poisson's equation with Dirichlet boundaries using the Thomas algorithm*/
void solvePotentialDirect(World &world)
{
	int ni = world.phi.size();	//number of mesh nodes
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
			d[i] = -world.rho[i]/EPS_0;
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
	world.phi[ni-1] = d[ni-1];
	for (int i=ni-2;i>=0;i--)
	{
		world.phi[i] = d[i] - c[i]*world.phi[i+1];
	}
}

/* solves potential using the Gauss Seidel Method, returns true if converged*/
bool solvePotentialGS(World &world, int max_it)
{
	double L2;
	double dx2 = world.dx*world.dx;		//precompute dx*dx
	const double w = 1.4;
	dvector &phi = world.phi;
	dvector &rho = world.rho;
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
				cout<<"Gauss-Seidel converged after "<<solver_it<<" iterations"<<endl;
				return false;
			}
		}
	}
	cout<<"Gauss-Seidel failed to converge, L2="<<L2<<endl;
	return true;
}
/* computes electric field by differentiating potential*/
void computeEF(World &world, bool second_order) {
	
	dvector &phi = world.phi;
	dvector &ef = world.ef;
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


bool Results::outputVTI() {

	ofstream out("results/field.vti");
	
	out<<"<VTKFile type=\"ImageData\">\n";
	out<<"<ImageData WholeExtent=\"0 "<<ni-1<<" 0 "<<nj-1<<" 0 0\" Origin=\"0 0 0\" Spacing=\"1 "<<1.5*ni/(double)nj<<" 1\">\n";
	out<<"<Piece Extent=\"0 "<<ni-1<<" 0 "<<nj-1<<" 0 0\">\n";
	out<<"<PointData>\n";

	out<<"<DataArray Name=\"phi\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	for (int u=0;u<ni*nj;u++) out<<phi[u]<<" ";
	out<<"\n</DataArray>\n";

	out<<"<DataArray Name=\"rho\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	for (int u=0;u<ni*nj;u++) out<<rho[u]<<" ";
	out<<"\n</DataArray>\n";

	out<<"<DataArray Name=\"ef\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	for (int u=0;u<ni*nj;u++) out<<ef[u]<<" ";
	out<<"\n</DataArray>\n";

	out<<"<DataArray Name=\"ndi\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	for (int u=0;u<ni*nj;u++) out<<ndi[u]<<" ";
	out<<"\n</DataArray>\n";

	out<<"<DataArray Name=\"nde\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	for (int u=0;u<ni*nj;u++) out<<nde[u]<<" ";
	out<<"\n</DataArray>\n";

	out<<"<DataArray Name=\"vi\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	for (int u=0;u<ni*nj;u++) out<<vi[u]<<" ";
	out<<"\n</DataArray>\n";

	out<<"<DataArray Name=\"ve\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	for (int u=0;u<ni*nj;u++) out<<ve[u]<<" ";
	out<<"\n</DataArray>\n";

	out<<"</PointData>\n";
	out<<"</Piece>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
	
	return true;
}

void Results::addData(Species &ions, Species &eles) {
	if (j>=nj) return;  // make sure we don't overflow the buffer 
	ions.computeVelocity();
	eles.computeVelocity();
	for (int i=0;i<ni;i++) this->phi[j*ni+i] = world.phi[i];
	for (int i=0;i<ni;i++) this->rho[j*ni+i] = world.rho[i];
	for (int i=0;i<ni;i++) this->ef[j*ni+i] = world.ef[i];
	for (int i=0;i<ni;i++) this->ndi[j*ni+i] = ions.den[i];
	for (int i=0;i<ni;i++) this->nde[j*ni+i] = eles.den[i];
	for (int i=0;i<ni;i++) this->vi[j*ni+i] = ions.vel[i];
	for (int i=0;i<ni;i++) this->ve[j*ni+i] = eles.vel[i];
	
	j++;
	
}
