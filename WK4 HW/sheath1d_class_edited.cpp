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
	//constructor
	World (int ni, double x0, double xm) : phi(ni), rho(ni), ef(ni) {
		this->ni = ni;
		this->x0 = x0;
		this->xm = xm;
		this->dx = (xm-x0)/(ni-1);
	}

	dvector phi;
	dvector rho;
	dvector ef;

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

double gather(double li, const dvector &f);
void scatter(double li, dvector &f, double val);

class Species {
public:
	// constructor
	//Species() {}
	//not sure why I needed to use () instead of {} for the den vector
	Species(World &world, double m, double q, double mpw, size_t N) : world{world}, den(world.ni) {
		(*this).m = m;
		 this->q = q;
		 this->mpw = mpw;
		 np = 0;
		 np_alloc = N;
		 part = new Particle[np_alloc];
	}
	dvector den;

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
  //compute average kinetic energy
  
double getAveKE() {
    double ave_KE;
    double sum;

        for (int i; i < np ; i++) {
            sum += part[i].v * part[i].v;
        }

    ave_KE = (0.5 * m * (sum / np) ) / 1.602E-19; //average particle KE in eV
    return ave_KE;
  }
  
  //Wall current
double getCurrent() {
	double time = world.dt * (world.ts - ts_last);
	double dnp = (np_last - np);

	double current;
	if (np_last >0) current = q * dnp * mpw / time; else current = 0;

	np_last = np;
	ts_last = world.ts;

	return current;
}

void computeNumberDensity() {
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

void rewindVelocity() {
	for (int p=0;p<np;p++) {
		dvector &ef = world.ef;
		double li = world.XtoL(part[p].x);       // get particle's grid index
		double E_p = gather(li, ef);            // electric field at particle position

		// here we go from -0.5 to +0.5
		part[p].v -= 0.5*q*E_p/m*world.dt;          // integrate velocity through dt
	}
}

void advanceSpecies() {
	dvector &ef = world.ef;

	for (int p=0;p<np;p++) {
		double li = world.XtoL(part[p].x);       // get particle's grid index

		double E_p = gather(li, ef);            // electric field at particle position

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

protected:
  int np_last = -1;
  double ts_last = -1;
  World &world;

public:

size_t np_alloc;
size_t np;
Particle *part;

};

class Results {
	public: 
	//constructor
	Results(World &world, int nj) : ni(world.ni), nj(nj), world{world} {
		
		
		//allocate bufferes and clear all bytes to 0
		
		phi = new double[world.ni*nj]; memset(phi,0,world.ni*nj*sizeof(double));
		rho = new double[world.ni*nj]; memset(rho,0,world.ni*nj*sizeof(double));
		ndi = new double[world.ni*nj]; memset(ndi,0,world.ni*nj*sizeof(double));
		nde = new double[world.ni*nj]; memset(nde,0,world.ni*nj*sizeof(double));
		
	}
		
	//destructor 
	~Results() {
		delete[] phi; phi = nullptr;
		delete[] rho; rho = nullptr;
		delete[] ndi; ndi = nullptr;
		delete[] nde; nde = nullptr;
		
	}

	//function decleration 
	//bool outputVTI();

	//databuffers
	double *phi = nullptr, *rho = nullptr;
	double *ndi = nullptr, *nde = nullptr;

	//dimensions 
	int ni = 0, nj = 0;

	bool outputVTI(Species &ions, Species &eles){
		dvector &phi = world.phi;
		dvector &rho = world.rho;
		dvector &ef = world.ef;
		dvector &ndi = ions.den;
		dvector &nde = eles.den;
		
		int width;
		int height;
		int depth;
		width = ni; 
		height = nj; //j axis coresponds to time
		depth = 0; //  = 0  with 2D data
		//spacing
		double dx, dy, dz;
		dx = 0.01;
		dy = 0.01;
		dz = 0;

		ofstream out("field.vti");

		 if (!out.is_open()) {
        	cerr << "Error: Could not open file for writing: field.vti" << endl;
        	return false;
    	}

		out << "<VTKFile type=\"ImageData\"> \n";
		out << "<ImageData WholeExtent = \"0 " << width -1 << " 0 " << height << " 0 " << depth  << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dy << " " << dz << "\">" << "\n";
		out << "<Piece Extent=\"0 " << width -1 << " 0 " << height-1 << " 0 " << depth  << "\">" << "\n";
		//out << "<Piece Extent=\"0 " << 20 << " 0 " << 20 << " 0 " << 14  << "\">" << "\n";
		out << "<PointData> \n";
		
		//phi
		out << "<DataArray type=\"Float64\" Name=\"phi\" NumberOfComponents=\"1\" format=\"ascii\">" << "\n";
		for (int u = 0; u < width *height; u++){
			out << phi[u] << " ";
			//out << "\n";
		}
		out << " \n </DataArray> \n";
		
		//rho
		out << "<DataArray type=\"Float64\" Name=\"rho\" NumberOfComponents=\"1\" format=\"ascii\">" << "\n";
		for (int u = 0; u < width *height; u++){
			out << rho[u] << " ";
			//out << "\n";
		}
		out << " \n </DataArray> \n";

		//nde
		out << "<DataArray type=\"Float64\" Name=\"nde\" NumberOfComponents=\"1\" format=\"ascii\">" << "\n";
		for (int u = 0; u < width *height; u++){
			out << nde[u] << " ";
			//out << "\n";
		}
		out << " \n </DataArray> \n";

		//ndi
		out << "<DataArray type=\"Float64\" Name=\"ndi\" NumberOfComponents=\"1\" format=\"ascii\">" << "\n";
		for (int u = 0; u < width *height; u++){
			out << ndi[u] << " ";
			//out << "\n";
		}
		out << " \n </DataArray> \n";


		//closing tags
		out << "</PointData> \n";
		out << "</Piece> \n";
		out << "</ImageData> \n";
		out << "</VTKFile> \n";

		out.close();
		return true;
	}

	//is this right? for step g, cleanup 2?
	void addData(Species &ions, Species &eles) {
		
		//do I need to also add if < 0 here?
		if ( j >= nj) {
			cerr << "Out of bounds" << endl;
			return;
		}

		int index = j * ni;
		for (int w = 0; w < ni; w++) {
			this->phi[index + w] = world.phi[w];
            this->rho[index + w] = world.rho[w];
            this->ndi[index + w] = ions.den[w];
            this->nde[index + w] = eles.den[w];
		}

		++j;
	}
    
	protected:
	World &world; //reference to world
	int j =0;


};

/*function prototypes*/
bool outputCSV(Species &ions, Species &eles, World &world);
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
	//Clean Up 1
	//World world;
	//world.ni = 41;	//number of nodes
	//world.x0 = 0;	//origin
	//world.xm = 0.1;	//opposite end
	//world.dx = (world.xm-world.x0)/(world.ni-1);	//node spacing
	//vector<double> phi(world.ni);
	//dvector rho(world.ni,QE*1e12);
	//dvector ef(world.ni);
	World world(41,0,0.1);

	//added in HW4 2 d)Hooks
	int num_ts = 20000;
	int add_skip = 100;
	Results results(world, num_ts/add_skip); //nj is the number of time entries

	constexpr int N = 100000;

	constexpr double n0 = 1e15;    // the desired number density in ... particles/m3 ?

	double n_real = n0*(world.xm-world.x0);   // number of real particles in domain of volume (xd-x0)*1*1

	Species ions = {world, 16*Const::AMU, Const::QE, n_real/N, N};
	Species eles(world, Const::ME, -Const::QE, n_real/N, N);

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

//do I need to delete these too for the
	vector<double> ndi(world.ni);
	dvector nde(world.ni);

	/*rewind for leapfrog!*/
		
	ions.advanceSpecies();
	eles.advanceSpecies();
	ions.computeNumberDensity();
	eles.computeNumberDensity();

	for (int i=0;i<world.ni;i++) {
		world.rho[i] = ndi[i]*ions.q + nde[i]*eles.q;
	}
	
	solvePotentialDirect(world);
	computeEF(world, true);
	ions.rewindVelocity();
	eles.rewindVelocity();

	
	// simulation main loop
    ofstream out("diag.csv");
	out << "something," << ions.getCurrent() << ", something else" << -eles.getCurrent() << endl;
    out << "Avg KE Electrons," << eles.getAveKE() << ", Avg KE Ions," << ions.getAveKE() << endl;
	out << "ts, time, num_ions, ion Current, num_electrons, elec Current" << endl;
	

	for (world.ts=0;world.ts< num_ts ;world.ts++) {

		ions.advanceSpecies();
		eles.advanceSpecies();

		ions.computeNumberDensity();
		eles.computeNumberDensity();

		for (int i=0;i<world.ni;i++) {
			world.rho[i] = ndi[i]*ions.q + nde[i]*eles.q;
		}

		solvePotentialDirect(world);
		computeEF(world, true);


		if (world.ts%100==0) {
            double sim_time;
            sim_time = world.ts * world.dt;
			out <<world.ts << "," << sim_time << "," << ions.np << "," << ions.getCurrent() << "," << eles.np << "," << -eles.getCurrent() << "\n";
			results.addData(ions,eles);
		}
		
		if (world.ts%100==0)
			//ouput to a CSV file for plotting
			outputCSV(ions,eles,world);

            
	}

	results.outputVTI(ions,eles);

	//cout<<"Done!"<<endl;


	return 0;	//normal exit
}

/*outputs the given fields to a CSV file, returns true if ok*/

//is this how I bring ions and eles in??
bool outputCSV(Species &ions, Species &eles, World &world)
{
    dvector &phi = world.phi;
	dvector &rho = world.rho;
	dvector &ef = world.ef;
	dvector &ndi = ions.den;
	dvector &nde = eles.den;

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

/*solves Poisson's equation with Dirichlet boundaries using the Thomas algorithm*/
void solvePotentialDirect(World &world)
{
	dvector &phi = world.phi;
	dvector &rho = world.rho;

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
bool solvePotentialGS(World &world, int max_it)
{
	dvector &phi = world.phi;
	dvector &rho = world.rho;
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
				cout<<"Gauss-Seidel converged after "<<solver_it<<" iterations"<<endl;
				return false;
			}
		}
	}
	cout<<"Gauss-Seidel failed to converge, L2="<<L2<<endl;
	return true;
}
/* computes electric field by differentiating potential*/

void computeEF(World &world, bool second_order)
{
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
