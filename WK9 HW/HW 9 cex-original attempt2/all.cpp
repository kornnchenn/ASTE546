#include <vector>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <memory>
#include <fstream>
#include <sstream>
#include <ostream>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

////////////////////World.h/////////////////////
//class Species;


//some typedefs
using Field = Field_<double>;
using FieldI = Field_<int>;
using Field3 = Field_<double3>;
using dvector = std::vector<double>;

/*define constants*/
namespace Const
{
	const double EPS_0 = 8.85418782e-12;  	// C/(V*m), vacuum permittivity
	const double QE = 1.602176565e-19;		// C, electron charge
	const double AMU = 1.660538921e-27;		// kg, atomic mass unit
	const double ME = 9.10938215e-31;		// kg, electron mass
	const double K = 1.380648e-23;			// J/K, Boltzmann constant
	const double PI = 3.141592653;			// pi
	const double EvToK = QE/K;				// 1eV in K ~ 11604
}

/*object for sampling random numbers*/
class Rnd {
public:
	//constructor: set initial random seed and distribution limits
	Rnd(): mt_gen{std::random_device()()}, rnd_dist{0,1.0} {}
	double operator() () {return rnd_dist(mt_gen);}

protected:
	std::mt19937 mt_gen;	    //random number generator
	std::uniform_real_distribution<double> rnd_dist;  //uniform distribution
};

extern Rnd rnd;		//tell the compiler that an object of type Rnd called rnd is defined somewhere

/*defines the computational domain*/
class World
{
public:	
	/*constructor, allocates memory*/
	World(int ni, int nj, int nk);

	/*functions to set mesh origin and spacing*/
	void setExtents(const double3 x0, const double3 xm);
	
	double3 getX0() const {return double3(x0);} 
	double3 getXm() const {return double3(xm);}
	double3 getXc() const {return double3(xc);}
	double3 getDh() const {return double3(dh);}
	double getCellVolume() const {return dh[0]*dh[1]*dh[2];} //cell volume

	/*functions for accessing time information*/
	int getTs() const {return ts;}
	double getTime() const {return time;}
	double getWallTime();  /*returns wall time in seconds*/
	double getDt() const {return dt;}
	bool isLastTimeStep() const {return ts==num_ts-1;}

	bool inBounds(double3 pos) {
		for (int i=0;i<3;i++)
			if (pos[i]<x0[i] || pos[i]>=xm[i]) return false;
		return true;
	}

	/*sets time step and number of time steps*/
	void setTime(double dt, int num_ts) {this->dt=dt;this->num_ts=num_ts;}
	
	/*advances to the next time step, returns true as long as more time steps remain*/
	bool advanceTime() {time+=dt;ts++;return ts<=num_ts;}

	/*checks and sets a steady state flag*/
	bool steadyState(std::vector<Species> &species);

	/*returns steady state flag*/
	bool isSteadyState() {return steady_state;}

	/*converts physical position to logical coordinate*/
	double3 XtoL(const double3 &x) const {
	  	double3 lc;
		lc[0] = (x[0]-x0[0])/dh[0];
		lc[1] = (x[1]-x0[1])/dh[1];
		lc[2] = (x[2]-x0[2])/dh[2];
		return lc;
	}

	/*returns cell index for a given point*/
	int3 XtoIJK(const double3 &x) const {
		double3 lc = XtoL(x);
		int3 ijk {(int)lc[0],(int)lc[1],(int)lc[2]};
		return ijk;
		}

	/*returns cell index for a position x*/
	int XtoC(double3 x) const {
		int3 i = XtoIJK(x);
		return i[2]*(nj-1)*(ni-1)+i[1]*(ni-1)+i[0];
	}


	/*converts logical coordinate to physical position*/
	double3 pos(const double3 &lc)
	{
		double3 x = x0+dh*lc;
		return x;
	}

	/*another form that takes 3 ints as inputs*/
	double3 pos(int i, int j, int k) {
		double3 x{(double)i,(double)j,(double)k};
		return pos(x);
	}

	int U(int i,int j, int k) {return object_id.U(i,j,k);}

	/*computes charge density from rho = sum(charge*den)*/
	void computeChargeDensity(std::vector<Species> &species);

	/*returns the system potential energy*/
	double getPE();

	/*sugarcubes a sphere centered at (x0,y0,z0)*/
	void addSphere(const double3 &x0, double radius, double phi_sphere);
	void addBox(const double3 &x0, const double3 &xm, double phi_sphere);
	
	/*return true if point x is inside or on the sphere*/
	bool inSphere(const double3 &x);
	bool inBox(const double3 &x);

	/*marks the k=0 face as Dirichlet*/
	//updated to make it a circle on the box max point (0.1 in this case)
	void addInlet(const double3 center, const double radius);

	//mesh geometry
	const int ni,nj,nk;	//number of nodes
	const int3 nn;	//another way to access node counts

	Field phi;			//potential
	Field rho;			//charge density
	Field node_vol;		//node volumes
	Field3 ef;			//electric field components
	FieldI object_id; 	//object id flag to flag fixed nodes

protected:
	double3 x0;	//mesh origin
	double3 dh;	//cell spacing

	double3 xm;		//origin-diagonally opposite corner (max bound)
	double3 xc;		//domain centroid

	double3 box_x0; // x0 point of rectangular object
	double3 box_xm; // xm point of rectangular object
	
	double dt = 0;		//time step length
	double time = 0;	//physical time
	int ts = -1;		//current time step
	int num_ts = 0;		//number of time steps

	/*sphere data*/
	double3 sphere_x0 {0,0,0};	//sphere centroid
	double sphere_rad2 = 0;		//sphere radius squared
	double3 _box_x0 {0,0,0}; //box corner
	double3 _box_xm {0,0,0}; //box corner

	std::chrono::time_point<std::chrono::high_resolution_clock> time_start;	//time at simulation start

	bool steady_state = false;	//set to true once steady state is reached
	double last_mass = 0;	//mass at the prior time step
	double last_mom = 0;	//momentum at the prior time step
	double last_en = 0;		//energy at the prior time step

	void computeNodeVolumes();
};

//////////////////////Field.h//////////////////////

template <typename T>
struct vec3 {
	vec3 (const T u, const T v, const T w) : d{u,v,w} {}
	vec3 (const T a[3]) : d{a[0],a[1],a[2]} {}
	vec3 (): d{0,0,0} {}
	T& operator[](int i) {return d[i];}
	T operator[](int i) const {return d[i];}
	vec3<T>& operator=(double s) {d[0]=s;d[1]=s;d[2]=s;return (*this);}
	vec3<T>& operator+=(vec3<T> o) {d[0]+=o[0];d[1]+=o[1];d[2]+=o[2];return(*this);}
	vec3<T>& operator-=(vec3<T> o) {d[0]-=o[0];d[1]-=o[1];d[2]-=o[2];return(*this);}
	vec3<T> operator/(double s) {vec3<T>o; o[0]=d[0]/s;o[1]=d[1]/s;o[2]=d[2]/s;return o;}
	vec3<T> operator/=(double s) {d[0]/=s;d[1]/=s;d[2]/=s;return (*this);}

	//dot product of two vectors
	friend T dot(const vec3<T> &v1, const vec3<T> &v2) {
		T s=0;	for (int i=0;i<3;i++) s+=v1[i]*v2[i];
		return s;	}

	//vector magnitude
	friend T mag(const vec3<T> &v) {return sqrt(dot(v,v));}

	//unit vector
	friend vec3<T> unit(const vec3<T> &v) {return vec3(v)/mag(v);}

	//cross product
	friend vec3<T> cross(const vec3<T> &a, const vec3<T> &b) {
		return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
	}

protected:
	T d[3];
};

//vec3-vec3 operations
template<typename T>	//addition of two vec3s
vec3<T> operator+(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a[0]+b[0],a[1]+b[1],a[2]+b[2]);	}
template<typename T>	//subtraction of two vec3s
vec3<T> operator-(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a[0]-b[0],a[1]-b[1],a[2]-b[2]);	}
template<typename T>	//element-wise multiplication of two vec3s
vec3<T> operator*(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a[0]*b[0],a[1]*b[1],a[2]*b[2]);	}
template<typename T>	//element wise division of two vec3s
vec3<T> operator/(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a[0]/b[0],a[1]/b[1],a[2]/b[2]);	}

//vec3 - scalar operations
template<typename T>		//scalar multiplication
vec3<T> operator*(const vec3<T> &a, T s) {
	return vec3<T>(a[0]*s, a[1]*s, a[2]*s);}
template<typename T>		//scalar multiplication 2
vec3<T> operator*(T s,const vec3<T> &a) {
	return vec3<T>(a[0]*s, a[1]*s, a[2]*s);}

//output
template<typename T>	//ostream output
std::ostream& operator<<(std::ostream &out, vec3<T>& v) {
	out<<v[0]<<" "<<v[1]<<" "<<v[2];
	return out;
}

using double3 = vec3<double>;
using int3 = vec3<int>;

template <typename T>
class Field_
{
public:
	
	/*constructor*/
	Field_(int ni, int nj, int nk) :
	ni{ni}, nj{nj}, nk{nk}
	{
		//allocate memory for a 3D array
		data = new T**[ni];
		for (int i=0;i<ni;i++)
		{
			data[i] = new T*[nj];
			for (int j=0;j<nj;j++) data[i][j] = new T[nk];
		}		

		clear();
	}

	//another constructor taking an int3
	Field_(int3 nn) : Field_(nn[0],nn[1],nn[2]) {};

	//copy constructor
	Field_(const Field_ &other):
	Field_{other.ni,other.nj,other.nk} {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++)
					data[i][j][k] = other(i,j,k);
	}
	
	//move constructor
	Field_(Field_ &&other):
		ni{other.ni},nj{other.nj},nk{other.nk} {
			data = other.data;	//steal the data
			other.data = nullptr;	//invalidate
	}

	//move assignment operator
	Field_& operator = (Field_ &&f) {data=f.data;
				f.data=nullptr; return *this;}

	//destructor: release memory
	~Field_() {
		//don't do anything if data is not allocated (or was moved away)
		if (data==nullptr) return;

		for (int i=0;i<ni;i++)
		{
			for (int j=0;j<nj;j++)
				delete[] data[i][j];

			delete[] data[i];
		}

		delete[] data;
	}

	//overloaded operator [] to allow direct access to data
	T** operator[] (int i) {return data[i];}

	/*returns data[i][j][k] marked as const to signal no data change*/
	T operator() (int i, int j, int k) const {return data[i][j][k];}

	/*sets all values to some scalar*/
	void operator =(double s) {
		for (int i=0;i<ni;i++)
		  for (int j=0;j<nj;j++)
		   for (int k=0;k<nk;k++)
			data[i][j][k] = s;
	  }

	/*performs element by element division by another field*/
	void operator /= (const Field_ &other) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++) {
					if (other.data[i][j][k]!=0)
					  data[i][j][k] /= other.data[i][j][k];
				else
					  data[i][j][k] = 0;
			  }
	}

	/*increments values by data from another field*/
	Field_& operator += (const Field_ &other) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++)
					data[i][j][k]+=other(i,j,k);
		return (*this);
	}

	/*performs element by element division by another field*/
	Field_& operator *= (double s) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++)
					data[i][j][k]*=s;
		return (*this);
	}

	//multiplication operator, returns f*s
	friend Field_<T> operator*(double s, const Field_<T>&f) {
		Field_<T> r(f);
		return r*=s;
	}

	//division of a field by a field of doubles
	friend Field_<T> operator*(const Field_<T>&f1, const Field_<T>&f2) {
		Field_<T> r(f1);
		for (int i=0;i<f1.ni;i++)
			for (int j=0;j<f1.nj;j++)
				for (int k=0;k<f1.nk;k++)
					r[i][j][k] = f1(i,j,k)*f2(i,j,k);
		return r;
	}

	//division of a field by a field of doubles
	friend Field_<T> operator/(const Field_<T>&f, const Field_<double>&d) {
		Field_<T> r(f);
		for (int i=0;i<f.ni;i++)
			for (int j=0;j<f.nj;j++)
				for (int k=0;k<f.nk;k++)
				{
					if (d(i,j,k)!=0)	//check for div by zero
						r[i][j][k] = f(i,j,k)/d(i,j,k);
					else
						r[i][j][k] = 0;
				}
		return r;
	}

	/*returns index for node (i,j,k)*/
	int U(int i, int j, int k) {return k*ni*nj+j*ni+i;}

	/*sets all data to zero*/
	void clear() {(*this)=0;}
	
	/* scatters scalar value onto a field at logical coordinate lc*/
	void scatter(double3 lc, T value)
	{
		int i = (int)lc[0];
		double di = lc[0]-i;
				
		int j = (int)lc[1];
		double dj = lc[1]-j;
		
		int k = (int)lc[2];
		double dk = lc[2]-k;
		
		data[i][j][k] += (T)value*(1-di)*(1-dj)*(1-dk);
		data[i+1][j][k] += (T)value*(di)*(1-dj)*(1-dk);
		data[i+1][j+1][k] += (T)value*(di)*(dj)*(1-dk);
		data[i][j+1][k] += (T)value*(1-di)*(dj)*(1-dk);
		data[i][j][k+1] += (T)value*(1-di)*(1-dj)*(dk);
		data[i+1][j][k+1] += (T)value*(di)*(1-dj)*(dk);
		data[i+1][j+1][k+1] += (T)value*(di)*(dj)*(dk);
		data[i][j+1][k+1] += (T)value*(1-di)*(dj)*(dk);		
	}

	/* gathers field value at logical coordinate lc*/
	T gather(double3 lc)
	{
		int i = (int)lc[0];
		double di = lc[0]-i;
				
		int j = (int)lc[1];
		double dj = lc[1]-j;
		
		int k = (int)lc[2];
		double dk = lc[2]-k;
					
		/*gather electric field onto particle position*/
		T val = data[i][j][k]*(1-di)*(1-dj)*(1-dk)+
				data[i+1][j][k]*(di)*(1-dj)*(1-dk)+
				data[i+1][j+1][k]*(di)*(dj)*(1-dk)+
				data[i][j+1][k]*(1-di)*(dj)*(1-dk)+
				data[i][j][k+1]*(1-di)*(1-dj)*(dk)+
				data[i+1][j][k+1]*(di)*(1-dj)*(dk)+
				data[i+1][j+1][k+1]*(di)*(dj)*(dk)+
				data[i][j+1][k+1]*(1-di)*(dj)*(dk);
				
		return val;
	}

	//incorporates new instantaneous values into a running average
	void updateAverage(const Field_ &I) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++)
					data[i][j][k] = (I(i,j,k)+ave_samples*data[i][j][k])/(ave_samples+1);
		++ave_samples;	//increment number of samples
	}

	template<typename S>
	friend std::ostream& operator<<(std::ostream &out, Field_<S> &f);
	const int ni,nj,nk;	//allocated dimensions

protected:
	T ***data;	/*data held by this field*/
	int ave_samples = 0;	//number of samples used for averaging
};

/*writes out data to a file stream*/
template<typename T>
std::ostream& operator<<(std::ostream &out, Field_<T> &f)
{
	for (int k=0;k<f.nk;k++,out<<"\n")
		for (int j=0;j<f.nj;j++)
			for (int i=0;i<f.ni;i++) out<<f.data[i][j][k]<<" ";
	return out;
}


/////////////////////Source.h//////////////////////
/*base class*/
class Source {
public:
	virtual void sample() = 0;
	virtual ~Source() {};	//to support destruction through base class
};


//simple monoenergetic source
class ColdBeamSource: public Source {
public:
	ColdBeamSource(Species &species, World &world, double v_drift, double den) :
		sp{species}, world{world}, v_drift{v_drift}, den{den} {}

	//generates particles
	void sample();

protected:
	Species &sp;	//reference to the injected species
	World &world;		//reference to world
	double v_drift;		//mean drift velocity
	double den;			//injection density
};

//simple monoenergetic source
class WarmBeamSource: public Source {
public:
	WarmBeamSource(Species &species, World &world, const double3 &x0, double rad, double mdot, double v_drift, double T) :
		sp{species}, world{world}, x0{x0}, rad{rad}, mdot{mdot}, v_drift{v_drift}, T{T} {}

	//generates particles
	void sample();

protected:
	Species &sp;	//reference to the injected species
	World &world;		//reference to world
	double v_drift;		//mean drift velocity
	double den;			//injection density
	double T;			//temperature
	double rad;			//radius of circle
	double mdot;		//mass flow rate of thruster
	double3 x0;         //center of source  (cirlce)
};

////////////////////Collisions.h//////////////////////

class Sigma {
public:
	virtual double operator() (double g) {return 0;};
	virtual ~Sigma () {};
};

class SigmaPoly: public Sigma {
	SigmaPoly (const std::vector<double> &coeffs) : coeffs{coeffs} {}

	double operator() (double g) {
		double r = coeffs[0];
		for (size_t i=1;i<coeffs.size();i++) {
			r+=coeffs[i]*g; g*=g;
		}
		return r;
	}


protected:
		std::vector<double> coeffs;
};


class Interaction {
public:
	virtual void apply(double dt) = 0;
	virtual ~Interaction() {}
};

//MCC Charge Exchange Collision
class MCC_CEX: public Interaction {
public:
	MCC_CEX(Species &source, Species &target, World &world) :
		source{source}, target{target}, world{world} {}
	void apply(double dt);
protected:
	Species &source;
	Species &target;
	World &world;
};

class DSMC_MEX: public Interaction {
public:
	DSMC_MEX(Species &species, World &world) : species{species}, world{world} {

		mr = species.mass*species.mass/(species.mass + species.mass);
		c[0] = 4.07e-10;
		c[1] = 0.77;
		c[2]= 2*Const::K*273.15/mr;	//Bird's reference params at 273.15 K
		c[3] = std::tgamma(2.5-c[1]); //Gamma(5/2-w)
	}
	void apply(double dt);
protected:
	void collide(double3 &vel1, double3 &vel2, double mass1, double mass2);
	double evalSigma(double g_rel) {
		return Const::PI*c[0]*c[0]*pow(c[2]/(g_rel*g_rel),c[1]-0.5)/c[3];

	}

	double sigma_cr_max = 1e-14;	//some initial value
	double mr;
	double c[4];
	Species &species;
	World &world;
};

class ChemistryIonize: public Interaction {
public:
	ChemistryIonize(Species &neutrals, Species &ions, World &world, double rate):
	neutrals{neutrals},ions{ions},world{world},rate{rate} { }
	void apply(double dt);
protected:
	Species &neutrals;
	Species &ions;
	World &world;
	double rate;
};

////////////////////Fields.h//////////////////////
template <typename T>
struct vec3 {
	vec3 (const T u, const T v, const T w) : d{u,v,w} {}
	vec3 (const T a[3]) : d{a[0],a[1],a[2]} {}
	vec3 (): d{0,0,0} {}
	T& operator[](int i) {return d[i];}
	T operator[](int i) const {return d[i];}
	vec3<T>& operator=(double s) {d[0]=s;d[1]=s;d[2]=s;return (*this);}
	vec3<T>& operator+=(vec3<T> o) {d[0]+=o[0];d[1]+=o[1];d[2]+=o[2];return(*this);}
	vec3<T>& operator-=(vec3<T> o) {d[0]-=o[0];d[1]-=o[1];d[2]-=o[2];return(*this);}
	vec3<T> operator/(double s) {vec3<T>o; o[0]=d[0]/s;o[1]=d[1]/s;o[2]=d[2]/s;return o;}
	vec3<T> operator/=(double s) {d[0]/=s;d[1]/=s;d[2]/=s;return (*this);}

	//dot product of two vectors
	friend T dot(const vec3<T> &v1, const vec3<T> &v2) {
		T s=0;	for (int i=0;i<3;i++) s+=v1[i]*v2[i];
		return s;	}

	//vector magnitude
	friend T mag(const vec3<T> &v) {return sqrt(dot(v,v));}

	//unit vector
	friend vec3<T> unit(const vec3<T> &v) {return vec3(v)/mag(v);}

	//cross product
	friend vec3<T> cross(const vec3<T> &a, const vec3<T> &b) {
		return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
	}

protected:
	T d[3];
};

//vec3-vec3 operations
template<typename T>	//addition of two vec3s
vec3<T> operator+(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a[0]+b[0],a[1]+b[1],a[2]+b[2]);	}
template<typename T>	//subtraction of two vec3s
vec3<T> operator-(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a[0]-b[0],a[1]-b[1],a[2]-b[2]);	}
template<typename T>	//element-wise multiplication of two vec3s
vec3<T> operator*(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a[0]*b[0],a[1]*b[1],a[2]*b[2]);	}
template<typename T>	//element wise division of two vec3s
vec3<T> operator/(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a[0]/b[0],a[1]/b[1],a[2]/b[2]);	}

//vec3 - scalar operations
template<typename T>		//scalar multiplication
vec3<T> operator*(const vec3<T> &a, T s) {
	return vec3<T>(a[0]*s, a[1]*s, a[2]*s);}
template<typename T>		//scalar multiplication 2
vec3<T> operator*(T s,const vec3<T> &a) {
	return vec3<T>(a[0]*s, a[1]*s, a[2]*s);}

//output
template<typename T>	//ostream output
std::ostream& operator<<(std::ostream &out, vec3<T>& v) {
	out<<v[0]<<" "<<v[1]<<" "<<v[2];
	return out;
}

using double3 = vec3<double>;
using int3 = vec3<int>;

template <typename T>
class Field_
{
public:
	
	/*constructor*/
	Field_(int ni, int nj, int nk) :
	ni{ni}, nj{nj}, nk{nk}
	{
		//allocate memory for a 3D array
		data = new T**[ni];
		for (int i=0;i<ni;i++)
		{
			data[i] = new T*[nj];
			for (int j=0;j<nj;j++) data[i][j] = new T[nk];
		}		

		clear();
	}

	//another constructor taking an int3
	Field_(int3 nn) : Field_(nn[0],nn[1],nn[2]) {};

	//copy constructor
	Field_(const Field_ &other):
	Field_{other.ni,other.nj,other.nk} {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++)
					data[i][j][k] = other(i,j,k);
	}
	
	//move constructor
	Field_(Field_ &&other):
		ni{other.ni},nj{other.nj},nk{other.nk} {
			data = other.data;	//steal the data
			other.data = nullptr;	//invalidate
	}

	//move assignment operator
	Field_& operator = (Field_ &&f) {data=f.data;
				f.data=nullptr; return *this;}

	//destructor: release memory
	~Field_() {
		//don't do anything if data is not allocated (or was moved away)
		if (data==nullptr) return;

		for (int i=0;i<ni;i++)
		{
			for (int j=0;j<nj;j++)
				delete[] data[i][j];

			delete[] data[i];
		}

		delete[] data;
	}

	//overloaded operator [] to allow direct access to data
	T** operator[] (int i) {return data[i];}

	/*returns data[i][j][k] marked as const to signal no data change*/
	T operator() (int i, int j, int k) const {return data[i][j][k];}

	/*sets all values to some scalar*/
	void operator =(double s) {
		for (int i=0;i<ni;i++)
		  for (int j=0;j<nj;j++)
		   for (int k=0;k<nk;k++)
			data[i][j][k] = s;
	  }

	/*performs element by element division by another field*/
	void operator /= (const Field_ &other) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++) {
					if (other.data[i][j][k]!=0)
					  data[i][j][k] /= other.data[i][j][k];
				else
					  data[i][j][k] = 0;
			  }
	}

	/*increments values by data from another field*/
	Field_& operator += (const Field_ &other) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++)
					data[i][j][k]+=other(i,j,k);
		return (*this);
	}

	/*performs element by element division by another field*/
	Field_& operator *= (double s) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++)
					data[i][j][k]*=s;
		return (*this);
	}

	//multiplication operator, returns f*s
	friend Field_<T> operator*(double s, const Field_<T>&f) {
		Field_<T> r(f);
		return r*=s;
	}

	//division of a field by a field of doubles
	friend Field_<T> operator*(const Field_<T>&f1, const Field_<T>&f2) {
		Field_<T> r(f1);
		for (int i=0;i<f1.ni;i++)
			for (int j=0;j<f1.nj;j++)
				for (int k=0;k<f1.nk;k++)
					r[i][j][k] = f1(i,j,k)*f2(i,j,k);
		return r;
	}

	//division of a field by a field of doubles
	friend Field_<T> operator/(const Field_<T>&f, const Field_<double>&d) {
		Field_<T> r(f);
		for (int i=0;i<f.ni;i++)
			for (int j=0;j<f.nj;j++)
				for (int k=0;k<f.nk;k++)
				{
					if (d(i,j,k)!=0)	//check for div by zero
						r[i][j][k] = f(i,j,k)/d(i,j,k);
					else
						r[i][j][k] = 0;
				}
		return r;
	}

	/*returns index for node (i,j,k)*/
	int U(int i, int j, int k) {return k*ni*nj+j*ni+i;}

	/*sets all data to zero*/
	void clear() {(*this)=0;}
	
	/* scatters scalar value onto a field at logical coordinate lc*/
	void scatter(double3 lc, T value)
	{
		int i = (int)lc[0];
		double di = lc[0]-i;
				
		int j = (int)lc[1];
		double dj = lc[1]-j;
		
		int k = (int)lc[2];
		double dk = lc[2]-k;
		
		data[i][j][k] += (T)value*(1-di)*(1-dj)*(1-dk);
		data[i+1][j][k] += (T)value*(di)*(1-dj)*(1-dk);
		data[i+1][j+1][k] += (T)value*(di)*(dj)*(1-dk);
		data[i][j+1][k] += (T)value*(1-di)*(dj)*(1-dk);
		data[i][j][k+1] += (T)value*(1-di)*(1-dj)*(dk);
		data[i+1][j][k+1] += (T)value*(di)*(1-dj)*(dk);
		data[i+1][j+1][k+1] += (T)value*(di)*(dj)*(dk);
		data[i][j+1][k+1] += (T)value*(1-di)*(dj)*(dk);		
	}

	/* gathers field value at logical coordinate lc*/
	T gather(double3 lc)
	{
		int i = (int)lc[0];
		double di = lc[0]-i;
				
		int j = (int)lc[1];
		double dj = lc[1]-j;
		
		int k = (int)lc[2];
		double dk = lc[2]-k;
					
		/*gather electric field onto particle position*/
		T val = data[i][j][k]*(1-di)*(1-dj)*(1-dk)+
				data[i+1][j][k]*(di)*(1-dj)*(1-dk)+
				data[i+1][j+1][k]*(di)*(dj)*(1-dk)+
				data[i][j+1][k]*(1-di)*(dj)*(1-dk)+
				data[i][j][k+1]*(1-di)*(1-dj)*(dk)+
				data[i+1][j][k+1]*(di)*(1-dj)*(dk)+
				data[i+1][j+1][k+1]*(di)*(dj)*(dk)+
				data[i][j+1][k+1]*(1-di)*(dj)*(dk);
				
		return val;
	}

	//incorporates new instantaneous values into a running average
	void updateAverage(const Field_ &I) {
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
				for (int k=0;k<nk;k++)
					data[i][j][k] = (I(i,j,k)+ave_samples*data[i][j][k])/(ave_samples+1);
		++ave_samples;	//increment number of samples
	}

	template<typename S>
	friend std::ostream& operator<<(std::ostream &out, Field_<S> &f);
	const int ni,nj,nk;	//allocated dimensions

protected:
	T ***data;	/*data held by this field*/
	int ave_samples = 0;	//number of samples used for averaging
};

/*writes out data to a file stream*/
template<typename T>
std::ostream& operator<<(std::ostream &out, Field_<T> &f)
{
	for (int k=0;k<f.nk;k++,out<<"\n")
		for (int j=0;j<f.nj;j++)
			for (int i=0;i<f.ni;i++) out<<f.data[i][j][k]<<" ";
	return out;
};


////////////////////Output.h//////////////////////
namespace Output {
	void fields(World &world, std::vector<Species> &species);
	void screenOutput(World &world, std::vector<Species> &species);
	void diagOutput(World &world, std::vector<Species> &species);
	void particles(World &world, std::vector<Species> &species, int num_parts);
};

//////////////////////PotentialSolver.h//////////////////////
//structure to hold data for a single row
template <int S>
struct Row {
	Row() {for (int i=0;i<S;i++) {a[i]=0;col[i]=-1;}}
	void operator= (const Row &o) {for (int i=0;i<S;i++) {a[i] = o.a[i];col[i]=o.col[i];}}
	double a[S];		//coefficients
	int col[S];
};

/*matrix with up to seven non zero diagonals*/
class Matrix
{
public:
    Matrix(int nr):nu{nr} {rows=new Row<nvals>[nr];}
    Matrix(const Matrix &o):Matrix(o.nu) {
    	for (int r=0;r<nu;r++) rows[r] = o.rows[r];
    };	//copy constructor
    ~Matrix() {if (rows) delete[] rows;}
	dvector operator*(dvector &v);	//matrix-vector multiplication

	double& operator() (int r, int c); //reference to A[r,c] value in a full matrix
	void clearRow(int r) {rows[r]=Row<nvals>();} //reinitializes a row
	Matrix diagSubtract(dvector &P);	//subtracts a vector from the diagonal
	Matrix invDiagonal();		//returns a matrix containing inverse of our diagonal
	double multRow(int r, dvector &x);	//multiplies row r with vector x

	static constexpr int nvals = 7;	//maximum 7 non-zero values
	const int nu;			//number of rows (unknowns)

protected:
	Row<nvals> *rows;	//row data
};

enum SolverType {GS, PCG, QN};

class PotentialSolver
{
public:
	/*constructor*/
	PotentialSolver(World &world, SolverType type, int max_it, double tol):
		world(world), solver_type(type), A(world.ni*world.nj*world.nk),
		max_solver_it(max_it), tolerance(tol) {
			buildMatrix();
		}

	/*sets reference values*/
	void setReferenceValues(double phi0, double Te0, double n0) {
		this->phi0 = phi0;
		this->Te0 = Te0;
		this->n0 = n0;
	}

	/*computes electric field = -gradient(phi)*/
	void computeEF();

	/*builds the "A" matrix for linear potential solver*/
    void buildMatrix();

    //calls the appropriate potential solver
    bool solve()
    {
    	switch(solver_type)
    	{
    		case GS: return solveGS();
    		case PCG: return solveNRPCG();
    		case QN: return solveQN();
    		default: return false;
    	}
    }
protected:
	World &world;
    SolverType solver_type;
    Matrix A;				//system matrix for the linear equation

    enum NodeType {REG,NEUMANN,DIRICHLET};
    std::vector<NodeType> node_type;	//flag for different node types

	unsigned max_solver_it;	//maximum number of solver iterations
	double tolerance;		//solver tolerance
	double phi0 = 0;		//reference plasma potential
	double n0 = 1e12;		//reference electron density
	double Te0 = 1.5;		//reference electron temperature in eV

	/*computes potential using quasineutral boltzmann model*/
	bool solveQN();

	/*solves non-linear potential using Gauss-Seidel*/
	bool solveGS();

	/*linear PCG solver for Ax=b system*/
	bool solvePCGLinear(Matrix &A, dvector &x, dvector &b);

	/*linear GS solver for Ax=b system*/
	bool solveGSLinear(Matrix &A, dvector &x, dvector &b);

	/*Newton Raphson solver for a nonlinear system, uses PCG for the linear solve*/
	bool solveNRPCG();
};

////////////////////Species.h//////////////////////
/** Data structures for particle storage **/
struct Particle
{
	double3 pos;			/*position*/
	double3 vel;			/*velocity*/
	double dt;				/*time step to push the particle through*/
	double mpw;				/*macroparticle weight*/
	
	Particle(double3 x, double3 v, double dt, double mpw):
	pos{x}, vel{v}, dt{dt}, mpw{mpw} { }
};

/*species container*/
class Species 
{
public:
	Species(std::string name, double mass, double charge, double mpw0, World &world) :
		name(name), mass(mass), charge(charge), mpw0(mpw0),
		den(world.nn), T(world.nn),	vel(world.nn),
		den_ave(world.nn),mpc(world.ni-1,world.nj-1,world.nk-1),
		n_sum(world.nn),nv_sum(world.nn),
		nuu_sum(world.nn), nvv_sum(world.nn), nww_sum(world.nn),
		world(world) { 	}

	/*returns the number of simulation particles*/
	size_t getNp()	{return particles.size();}

	/*returns the number of real particles*/
	double getRealCount();

	/*returns the species momentum*/
	double3 getMomentum();

	/*returns the species kinetic energy*/
	double getKE();

	/*moves all particles using electric field ef[]*/
	void advance();

	/*compute number density*/
	void computeNumberDensity();

	/*samples velocity moments*/
	void sampleMoments();

	/*uses sampled data to compute velocity and temperature*/
	void computeGasProperties();

	/*computes number of macroparticles per cell*/
	void computeMPC();

	/*clears sampled moment data*/
	void clearSamples();

	/*adds a new particle*/
	void addParticle(double3 pos, double3 vel) {addParticle(pos,vel,world.getDt());}
	void addParticle(double3 pos, double3 vel, double dt) {addParticle(pos,vel,dt,mpw0);}
	void addParticle(double3 pos, double3 vel, double dt, double mpw);

	/*updates number density*/
	void updateAverages() {den_ave.updateAverage(den);}

	/*returns random thermal velocity*/
	double sampleVth(double T);

	/*samples random isotropic velocity*/
	double3 sampleIsotropicVel(double T);

	const std::string name;			/*species name*/
	const double mass;			/*particle mass in kg*/
	const double charge;		/*particle charge in Coulomb*/
	const double mpw0;			/*default macroparticle weight*/
	
	std::vector<Particle> particles;	/*contiguous array for storing particles*/
	Field den;			/*number density*/
	Field T;			/*temperature*/
	Field3 vel;			/*stream velocity*/
	Field den_ave;		/*averaged number density*/
	Field mpc;			/*macroparticles per cell*/

protected:

	Field n_sum;
	Field3 nv_sum;
	Field nuu_sum,nvv_sum,nww_sum;
	World &world;
};


////////////////////Main.cpp//////////////////////
using namespace std;		//to avoid having to write std::cout
using namespace Const;		//to avoid having to write Const::ME

/*program execution starts here*/
int main(int argc, char *args[])
{
	/*initialize domain*/
    World world(31,81,71);
    world.setExtents({0,-0.4,0},{0.3,0.4,0.7});
    world.setTime(2e-7,5000);


	/*set objects*/
	double phi_sphere = -100;		//set default
	if (argc>1)
		phi_sphere = atof(args[1]);	//convert argument to float
	cout<<"Sphere potential: "<<phi_sphere<<" V"<<endl;
    world.addBox({-0.1,-0.1,0},{0.1,0.1,0.1},phi_sphere);


	/*set up particle species*/
    vector<Species> species;
    species.push_back(Species("X", 131.3*AMU, 0, 2e9, world));
	species.push_back(Species("X+", 13*AMU, QE, 5e8, world));
	Species &neutrals = species[0];
	Species &ions = species[1];
	
	/*setup injection sources
	const double nda = 2e20;			//neutral density
	const double ndi = 1e10;			//mean ion density
	vector<unique_ptr<Source>> sources;
	sources.emplace_back(new WarmBeamSource(neutrals,world,7000,nda,1000));	*/	//neutral source

	/*setup injection sources*/
    //const double nda = 2e20;            //neutral density
    const double ndi = 1e17;            //mean ion density --- changed in step 5
    vector<unique_ptr<Source>> sources;
    //sources.emplace_back(new WarmBeamSource(neutrals,world,{0,0,0.1}, 0.02, 0.1e-6 ,500, 300));  //neutral source
    sources.emplace_back(new WarmBeamSource(ions,world,{0,0,0.1}, 0.02, 0.8e-6 , 15000, 1e4));   //ions source

	
	/*setup material interactions*/
	vector<unique_ptr<Interaction>> interactions;
	interactions.emplace_back(new MCC_CEX(ions,neutrals,world));

 /*initialize potential solver and solve initial potential*/
    PotentialSolver solver(world,SolverType::QN,1000,1e-4);
    solver.setReferenceValues(30,2.5,ndi); //30 V, 2.5eV
    solver.solve();


    /*obtain initial electric field*/
    solver.computeEF();

    /* main loop*/
	while(world.advanceTime())
    {
		/*inject particles*/
    	for (auto &source:sources)
    		source->sample();

    	/*perform material interactions*/
    	for (auto &interaction:interactions)  interaction->apply(world.getDt());

		/*move particles*/
		for (Species &sp:species)
		{
			sp.advance();
			sp.computeNumberDensity();
			sp.sampleMoments();
			sp.computeMPC();
		}

		/*compute charge density*/
		world.computeChargeDensity(species);

        /*update potential*/
       // solver.solve();

        /*obtain electric field*/
        solver.computeEF();

        /*update averages at steady state*/
        if (world.steadyState(species)) {
        	for (Species &sp:species)
        		sp.updateAverages();
        }

		/*screen and file output*/
        Output::screenOutput(world,species);
        Output::diagOutput(world,species);

		/*periodically write out results*/
        if (world.getTs()%50==0 || world.isLastTimeStep()) {
			Output::fields(world, species);
		//	Output::particles(world, species,10000);
        }

    }
	
	/* grab starting time*/
	cout<<"Simulation took "<<world.getWallTime()<<" seconds"<<endl;
	return 0;		//indicate normal exit
}


//////////////Collisions.cpp////////////////////
//approximates ionization, not mass conserving
void ChemistryIonize::apply(double dt) {
	double dV = world.getCellVolume();
	//loop over all cells
	for (int i=0;i<world.ni-1;i++)
		for (int j=0;j<world.nj-1;j++)
			for (int k=0;k<world.nk-1;k++) {
				//evaluate electron and neutral density at cell center
				double na = neutrals.den.gather({i+0.5,j+0.5,k+0.5});
				double ne = world.rho.gather({i+0.5,j+0.5,k+0.5})/Const::QE; //assume QN
				double dni = rate*na*ne*dt;

				/*number of macroparticles to create*/
				int num_p = (int)(dni*dV/ions.mpw0 + rnd());
				for (int p=0;p<num_p;p++) {
					/*sample a random particle in the cell*/
					double3 pos = world.pos(i+rnd(),j+rnd(),k+rnd());
					double3 vel {0,0,0}; //should sample from Maxwellian instead
					ions.addParticle(pos,vel,ions.mpw0);
				}
			

			}
}

//return 1e-18*pow(rel_g,-0.5);

/*contant cross-section for the demo*/
double evalSigma(double rel_g)
{
	return 1e-16;
}


/*MCC*/
void MCC_CEX::apply(double dt)
{
		/*set pointers to target data*/
		Field &target_den = target.den;
		Field &target_temp = target.T;
		Field3 &target_vel = target.vel;

		/*loop over all particles*/
        for (Particle &part:source.particles)
        {
			/*get target velocity*/
			double3 lc = world.XtoL(part.pos);
            			
			double3 vt = target_vel.gather(lc);

            /*get target density*/
            double nn = target_den.gather(lc);
           
			/*compute cross-section, function of relative velocity*/
            double3 v_rel = part.vel-vt;
			double v_rel_mag = mag(v_rel);
            
            /*evaluate cross-section */
            //double sigma = 1e-16;
			double sigma = 5e-18;
            
            /*compute probability*/
            double P = 1 - exp(-nn*sigma*v_rel_mag*dt);
            
            /*compare to a random number to see if collision happened*/
            if (P>=rnd())
            {
                /*sample a virtual target particle*/
                double T_target = target_temp.gather(lc);	//sample target temperature
                double3 v_target = target.sampleIsotropicVel(T_target);
                //part.vel = v_target; 	//CEX collision
                
                part.vel = 0;
            }
        }
    }
    

/* collides two particles*/
void DSMC_MEX::collide(double3 &vel1, double3 &vel2, double mass1, double mass2)
{
	double3 cm = (mass1*vel1 + mass2*vel2)/(mass1+mass2);

	/*relative velocity, magnitude remains constant through the collision*/
	double3 cr = vel1 - vel2;
	double cr_mag = mag(cr);

	/*pick two random angles, per Bird's VHS method*/
	double cos_chi = 2*rnd()-1;
	double sin_chi = sqrt(1-cos_chi*cos_chi);
	double eps = 2*Const::PI*rnd();

	/*perform rotation*/
	cr[0] = cr_mag*cos_chi;
	cr[1] = cr_mag*sin_chi*cos(eps);
	cr[2] = cr_mag*sin_chi*sin(eps);

	/*post collision velocities*/
	vel1 = cm + mass2/(mass1+mass2)*cr;
	vel2 = cm - mass1/(mass1+mass2)*cr;
}


/*DSMC*/
void DSMC_MEX::apply(double dt)
{
	/*first we need to sort particles to cells*/
	vector<Particle*> *parts_in_cell;
	int n_cells = (world.ni-1)*(world.nj-1)*(world.nk-1);
	parts_in_cell = new vector<Particle*> [n_cells];

	/*sort particles to cells*/
	for (Particle &part:species.particles)
	{
		int c = world.XtoC(part.pos);
		parts_in_cell[c].push_back(&part);
	}

	double sigma_cr_max_temp = 0;	/*reset for max computation*/
	double dV = world.getCellVolume();	/*internal cell volume*/
	double Fn = species.mpw0;	/*specific weight, using Bird's notation*/
	int num_cols=0;	/*reset collision counter*/
					
	/*now perform collisions*/
	for (int c=0;c<n_cells;c++)
	{
		vector<Particle*> &parts = parts_in_cell[c];
		int np = parts.size();
		if (np<2) continue;

		/*compute number of groups according to NTC*/
		double ng_f = 0.5*np*np*Fn*sigma_cr_max*dt/dV;
		int ng = (int)(ng_f+0.5);	/*number of groups, round*/
	
		/*assumes at least two particles per cell*/
		for (int g=0;g<ng;g++)
		{
			int p1, p2;
			p1 = (int)(rnd()*np);		/*returns some number between 0 and np-1 inclusive*/
		
			do {
				p2 = (int)(rnd()*np);
			} while (p2==p1);

			/*compute relative velocity*/
			double3 cr_vec = parts[p1]->vel - parts[p2]->vel;
			double cr = mag(cr_vec);

			/*evaluate cross section*/
			double sigma = evalSigma(cr);

			/*eval sigma_cr*/
			double sigma_cr=sigma*cr;

			/*update sigma_cr_max*/
			if (sigma_cr>sigma_cr_max_temp)
				sigma_cr_max_temp=sigma_cr;

			/*eval prob*/
			double P=sigma_cr/sigma_cr_max;

			/*did the collision occur?*/
			if (P>rnd())
			{
				num_cols++;
				collide(parts[p1]->vel,parts[p2]->vel,species.mass, species.mass);
			}
		}
	}

	delete[] parts_in_cell;

	if (num_cols){
		sigma_cr_max = sigma_cr_max_temp;
	}

}

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

////////////////////Species.cpp//////////////////////
/*updates velocities and positions of all particles of this species*/
void Species::advance()
{
		/*loop over all particles*/
	for (Particle &part: particles)
	{
		/*increment particle's dt by world dt*/
		part.dt += world.getDt();

		/*get logical coordinate of particle's position*/
		double3 lc = world.XtoL(part.pos);
		
		/*electric field at particle position*/
		double3 ef_part = world.ef.gather(lc);
			
		/*update velocity from F=qE*/
		part.vel += ef_part*(part.dt*charge/mass);

		/*update position from v=dx/dt, take into account particle bounces*/
		int n_bounces = 0;

		/*keep iterate while time remains and the particle is alive*/
		while (part.dt>0 && part.mpw>0) {
			double3 pos_old = part.pos;
			part.pos += part.vel*part.dt;

			/*did this particle leave the domain?*/
			if (!world.inBounds(part.pos) || world.inBox(part.pos))
			{
				part.mpw = 0;	//kill the particle
			}

			//this particle finished the whole step
			part.dt = 0;

			//kill stuck particles
			if (++n_bounces>20) {std::cerr<<"Stuck particle!"<<std::endl;part.mpw = 0;}
		}
	}

	/*perform a particle removal step, dead particles are replaced by the entry at the end*/
	size_t np = particles.size();
	for (size_t p=0;p<np;p++)
	{
		if (particles[p].mpw>0) continue;	//ignore live particles
		particles[p] = particles[np-1]; //fill the hole
		np--;	//reduce count of valid elements
		p--;	//decrement p so this position gets checked again
	}

	//now delete particles[np:end]
	particles.erase(particles.begin()+np,particles.end());
}

/*adds a new particle, rewinding velocity by half dt*/
void Species::addParticle(double3 pos, double3 vel, double dt, double mpw)
{
	//don't do anything (return) if pos outside domain bounds [x0,xd)
	if (!world.inBounds(pos)) return;

	//get particle logical coordinate
	double3 lc = world.XtoL(pos);
	
	//evaluate electric field at particle position
    double3 ef_part = world.ef.gather(lc);

	//rewind velocity back by 0.5*dt*ef
    vel -=  charge/mass*ef_part*(0.5*world.getDt());

    //add to list
    particles.emplace_back(pos,vel,dt,mpw);
}

/*returns the number of real particles*/
double Species::getRealCount() {
	double mpw_sum = 0;
	for (Particle &part:particles)
		mpw_sum+=part.mpw;
	return mpw_sum;
}

/* returns the species momentum*/
double3 Species::getMomentum() {
	double3 mom;
	for (Particle &part:particles)
		mom+=part.mpw*part.vel;
	return mass*mom;
}

/* returns the species kinetic energy*/
double Species::getKE() {
	double ke = 0;
	for (Particle &part:particles)
	{
		double v2 = part.vel[0]*part.vel[0] + part.vel[1]*part.vel[1] + part.vel[2]*part.vel[2];
		ke += part.mpw*v2;
	}
	return 0.5*mass*ke;
}

/*returns random thermal velocity*/
double Species::sampleVth(double T)
{
	//thermal velocity
	double v_th = sqrt(2*Const::K*T/mass);
	//get three random velocity components
	double v1 = v_th*(rnd()+rnd()+rnd()-1.5);
	double v2 = v_th*(rnd()+rnd()+rnd()-1.5);
	double v3 = v_th*(rnd()+rnd()+rnd()-1.5);
	return 3/sqrt(2+2+2)*sqrt(v1*v1+v2*v2+v3*v3);	//magnitude
}

/*returns random isotropic velocity*/
double3 Species::sampleIsotropicVel(double T) {
	double theta = 2*Const::PI*rnd();
    double r = -1.0+2*rnd();    //pick a random direction for d[0]
    double a = sqrt(1-r*r); //scaling for unity magnitude

    double3 d;
    d[0] = r;
	d[1] = cos(theta)*a;
    d[2] = sin(theta)*a;

    double v_th = sampleVth(T);
    double3 vel = v_th*d;
    return vel;
}

/*compute number density*/
void Species::computeNumberDensity()
{
	den.clear();
	for (Particle &part:particles)
	{
		double3 lc = world.XtoL(part.pos);
		den.scatter(lc, part.mpw);
	}

	//divide by node volume
	den /= world.node_vol;
}

/*samples velocity moments*/
void Species::sampleMoments() {
	for (Particle &part:particles)
	{
		double3 lc = world.XtoL(part.pos);
		n_sum.scatter(lc, part.mpw);
		nv_sum.scatter(lc,part.mpw*part.vel);
		nuu_sum.scatter(lc,part.mpw*part.vel[0]*part.vel[0]);
		nvv_sum.scatter(lc,part.mpw*part.vel[1]*part.vel[1]);
		nww_sum.scatter(lc,part.mpw*part.vel[2]*part.vel[2]);
	}
}

/*uses sampled data to compute velocity and temperature*/
void Species::computeGasProperties() {
	vel = nv_sum/n_sum;	//stream velocity

	for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
			for (int k=0;k<world.nk;k++) {
				double count = n_sum(i,j,k);
				if (count<=0) {T[i][j][k] = 0; continue;}

				double u_ave = vel(i,j,k)[0];
				double v_ave = vel(i,j,k)[1];
				double w_ave = vel(i,j,k)[2];
				double u2_ave = nuu_sum(i,j,k)/count;
				double v2_ave = nvv_sum(i,j,k)/count;
				double w2_ave = nww_sum(i,j,k)/count;

				double uu = u2_ave - u_ave*u_ave;
				double vv = v2_ave - v_ave*v_ave;
				double ww = w2_ave - w_ave*w_ave;
				T[i][j][k] = mass/(2*Const::K)*(uu+vv+ww);
			}
}

/*computes number of macroparticles per cell*/
void Species::computeMPC() {
	mpc.clear();
	for (Particle &part:particles) {
		int3 ijk = world.XtoIJK(part.pos);
		int i = ijk[0], j = ijk[1], k=ijk[2];
		mpc[i][j][k] += 1;
	}
}

/*clears sampled moment data*/
void Species::clearSamples() {
	n_sum = 0; nv_sum = 0; nuu_sum = 0; nvv_sum = 0; nww_sum=0;
}


//samples monoenergetic particles according to a prescribed density
void ColdBeamSource::sample()
{
	double3 dh = world.getDh();
	double3 x0 = world.getX0();

	//area of the XY plane, A=Lx*Ly
	double Lx = dh[0]*(world.ni-1);
	double Ly = dh[1]*(world.nj-1);
	double A = Lx*Ly;

	//compute number of real particles to generate: (#/s) = n*v*A; # = (#/s)*dt
	double num_real = den*v_drift*A*world.getDt();

	//number of simulation particles
	int num_sim = (int)(num_real/sp.mpw0+rnd());

	//inject particles
	for (int i=0;i<num_sim;i++)
	{
		double3 pos {x0[0]+rnd()*Lx, x0[1]+rnd()*Ly, x0[2]};
		double3 vel {0,0,v_drift};
		sp.addParticle(pos,vel);
	}
}

//samples particles with finite thermal and drift velocity
void WarmBeamSource::sample()
{
	double3 dh = world.getDh(); //cell spacing
	double3 x0 = world.getX0(); //mesh origin

	//area of the XY plane, A=Lx*Ly
	//double Lx = dh[0]*(world.ni-1);
	//double Ly = dh[1]*(world.nj-1);
	//double A = Lx*Ly;

	//compute number of real particles to generate: (#/s) = n*v*A; # = (#/s)*dt
	//double num_real = den*v_drift*A*world.getDt();
	
	//compute number of real particles to generate (#/s)
	//mass flow rate = mdot, divided by atomic mass (kg/particle) == particles/sec
	double num_real = mdot / (sp.mass); 

	//number of simulation particles
	int num_sim = (int)(num_real/sp.mpw0+rnd());

	//inject particles
	for (int i=0;i<num_sim;i++)
	{
		//double3 pos {x0[0]+rnd()*Lx, x0[1]+rnd()*Ly, x0[2]};
		double3 pos;
		double angle;
		double radius;
		double x, y;
		angle = rnd() * 360;
		radius = sqrt(rnd() * (0.02 * 0.02));
		x = radius * cos(angle*(M_PI/180));
		y = radius * sin(angle*(M_PI/180));
		//pos = {x,y, x0[2]};
		pos = {x,y, 0.1}; // if I use the z fromt the origin, x0[2] then won't that be at 0 and not coming out of the box?
		double3 vel = sp.sampleIsotropicVel(T);
        vel[2] += v_drift; //add drift component

        sp.addParticle(pos,vel);
	}
}

////////////////////PotentialSolver.cpp//////////////////////

//matrix-vector multiplication
dvector Matrix::operator*(dvector &v) {
	dvector r(nu);
	for (int u=0;u<nu;u++) {
		auto &row = rows[u];
		r[u] = 0;
		for (int i=0;i<nvals;i++){
			if (row.col[i]>=0) r[u]+=row.a[i]*v[row.col[i]];
			else break;	//end at the first -1
		}
	}
	return r;
}

//returns reference to A[r,c] element in the full matrix
double& Matrix::operator()(int r, int c){
	//find this entry
	auto &row = rows[r]; int v;
	for (v=0;v<nvals;v++)
	{
		if (row.col[v]==c) break;	//if found
		if (row.col[v]<0) {row.col[v]=c;   //set
						   break;}
	}
	assert(v!=nvals);	//check for overflow
	return row.a[v];
}

/*returns inverse of a diagonal preconditioner*/
Matrix Matrix::invDiagonal()
{
	Matrix M(nu);
	for (int r=0;r<nu;r++)	M(r,r) = 1.0/(*this)(r,r);
   return M;
}

/*subtracts diagonal matrix diag from A*/
Matrix Matrix::diagSubtract(dvector &P) {
	Matrix M(*this);	//make a copy
	for (int u=0;u<nu;u++) M(u,u)=(*this)(u,u)-P[u];
	return M;
}

//multiplies row r with vector x
double Matrix::multRow(int r, dvector &x){
	auto &row = rows[r];
	double sum=0;
	for (int i=0;i<nvals;i++)
	{
		if (row.col[i]>=0) sum+=row.a[i]*x[row.col[i]];
		else break;
	}
	return sum;
}


dvector operator-(const dvector &a, const dvector &b) {
	size_t nu = a.size();
	dvector r(nu);
	for (size_t u=0;u<nu;u++) r[u] = a[u]-b[u];
	return r;
}

dvector operator+(const dvector &a, const dvector &b) {
	size_t nu = a.size();
	dvector r(nu);
	for (size_t u=0;u<nu;u++) r[u] = a[u]+b[u];
	return r;
}

dvector operator*(const double s, const dvector &a) {
	size_t nu = a.size();
	dvector r(nu);
	for (size_t u=0;u<nu;u++) r[u] = s*a[u];
	return r;
}

/*vector math helper functions*/
namespace vec
{
	/*returns sum of v1[i]*v2[i]*/
	double dot(dvector v1, dvector v2)
	{
	    double dot = 0;
	    size_t nu = v1.size();
        for (size_t j=0;j<nu;j++)
            dot+=v1[j]*v2[j];
        return dot;
	}

	/*returns l2 norm*/
	double norm(dvector v)
	{
		double sum = 0;
		int nu = v.size();
        for (int j=0;j<nu;j++)
            sum+=v[j]*v[j];
		return sqrt(sum/nu);
	}

	/** converts 3D field to a 1D vector*/
	dvector deflate(Field &f3)
	{
		dvector r(f3.ni*f3.nj*f3.nk);
		for (int i=0;i<f3.ni;i++)
			  for (int j=0;j<f3.nj;j++)
					for (int k=0;k<f3.nk;k++)
						 r[f3.U(i,j,k)] = f3[i][j][k];
		return r;
	}

	/** converts 1D vector to 3D field*/
	void inflate(dvector &d1, Field& f3)
	{
		for (int i=0;i<f3.ni;i++)
			for (int j=0;j<f3.nj;j++)
				for (int k=0;k<f3.nk;k++)
					f3[i][j][k] = d1[f3.U(i,j,k)];
	}

};

//constructs the coefficient matrix
void PotentialSolver::buildMatrix()
{
	double3 dh = world.getDh();
	double idx = 1.0/dh[0];
	double idy = 1.0/dh[1];
	double idz = 1.0/dh[2];
    double idx2 = idx*idx;	/*1/(dx*dx)*/
	double idy2 = idy*idy;
	double idz2 = idz*idz;
	int ni = world.ni;
	int nj = world.nj;
	int nk = world.nk;
	int nu = ni*nj*nk;

	/*make sure we have space for node types*/
	node_type.resize(nu);

	/*solve potential*/
	for (int k=0;k<nk;k++)
        for (int j=0;j<nj;j++)
        	for (int i=0;i<ni;i++)
            {
                int u = world.U(i,j,k);
                A.clearRow(u);
                //dirichlet node?
				if (world.object_id[i][j][k]>0)
                {
                    A(u,u)=1;	//set 1 on the diagonal
                    node_type[u] = DIRICHLET;
                    continue;
                }

				//Neumann boundaries
				node_type[u] = NEUMANN;		//set default
                if (i==0) {A(u,u)=idx;A(u,u+1)=-idx;}
                else if (i==ni-1) {A(u,u)=idx;A(u,u-1)=-idx;}
                else if (j==0) {A(u,u)=idy;A(u,u+ni)=-idy;}
                else if (j==nj-1) {A(u,u)=idy;A(u,u-ni)=-idy;}
                else if (k==0) {A(u,u)=idz;A(u,u+ni*nj)=-idz;}
				else if (k==nk-1) {
					A(u,u)=idz;
					A(u,u-ni*nj)=-idz;}
                else {
                	//standard internal stencil
                	A(u,u-ni*nj) = idz2;
                	A(u,u-ni) = idy2;
                	A(u,u-1) = idx2;
                	A(u,u) = -2.0*(idx2+idy2+idz2);
                	A(u,u+1) = idx2;
                	A(u,u+ni) = idy2;
                	A(u,u+ni*nj) = idz2;
                	node_type[u] = REG;	//regular internal node
                }
            }
	solveQN();
}

/*quasi-neutral potential solver*/
bool PotentialSolver::solveQN()
{
	Field& phi = world.phi;
    Field& rhoi = world.rho;
    double rho0 = n0*QE;
    double rho_ratio_min = 1e-6;

    for (int i=0;i<world.ni;i++)
        for (int j=0;j<world.nj;j++)
            for (int k=0;k<world.nk;k++)
            {
                if (world.object_id[i][j][k]>0) continue; /*skip Dirichlet nodes*/

				double rho_ratio = rhoi[i][j][k]/rho0;
				if (rho_ratio<rho_ratio_min) rho_ratio=rho_ratio_min;
				phi[i][j][k] = phi0 + Te0*log(rho_ratio);
            }
    return true;
}

/*Newton Raphson solver for a nonlinear system, using PCG for the linear solve	*/
bool PotentialSolver::solveNRPCG()
{
	/*main NR iteration loop*/
	const int NR_MAX_IT=20;		/*maximum number of NR iterations*/
	const double NR_TOL = 1e-3;
	int nu = A.nu;

	Matrix J(nu);
	dvector P(nu);
	dvector y(nu);
	dvector x = vec::deflate(world.phi);
	dvector b = vec::deflate(world.rho);

	/*set RHS to zero on boundary nodes (zero electric field)
      and to existing potential on fixed nodes */
    for (int u=0;u<nu;u++)
    {
		if (node_type[u]==NEUMANN) b[u] = 0;			/*neumann boundary*/
        else if (node_type[u]==DIRICHLET) b[u] = x[u];	/*dirichlet boundary*/
        else b[u] = -b[u]/EPS_0;            /*regular node*/
    }

	double norm;
	bool converged=false;
	for(int it=0;it<NR_MAX_IT;it++)
	{
		/*compute F by first subtracting the linear term */
		dvector F = A*x-b;

		/*subtract b(x) on regular nodes*/
		for (int n=0;n<nu;n++)
			if (node_type[n]==REG)	/*regular nodes*/
				F[n] -= QE*n0*exp((x[n]-phi0)/Te0)/EPS_0;

		/*Compute P, diagonal of d(bx)/dphi*/
		for (int n=0;n<nu;n++)
		{
			if (node_type[n]==REG)
				P[n] = n0*QE/(EPS_0*Te0)*exp((x[n]-phi0)/Te0);
		}

		/*Compute J = A-diag(P)*/
		Matrix J = A.diagSubtract(P);

		/*solve Jy=F*/
		if (!solvePCGLinear(J,y,F))
			solveGSLinear(J,y,F);

		/*clear any numerical noise on Dirichlet nodes*/
		for (int u=0;u<nu;u++)
			if (node_type[u]==DIRICHLET) y[u]=0;

		/*x=x-y*/
		x = x-y;

		norm=vec::norm(y);
		//cout<<"NR norm: "<<norm<<endl;

		if (norm<NR_TOL)
		{
			converged=true;
			break;
		}
	}

	if (!converged)
		cout<<"NR+PCG failed to converge, norm = "<<norm<<endl;

	/*convert to 3d data*/
	vec::inflate(x,world.phi);
	return converged;
}

/*PCG solver for a linear system Ax=b*/
bool PotentialSolver::solvePCGLinear(Matrix &A, dvector &x, dvector &b)
{
	bool converged= false;

	double l2 = 0;
	Matrix M = A.invDiagonal(); //inverse of Jacobi preconditioner

	/*initialization*/
	dvector g = A*x-b;
	dvector s = M*g;
	dvector d = -1*s;

	for (unsigned it=0;it<max_solver_it;it++)
	{
		dvector z = A*d;
		double alpha = vec::dot(g,s);
		double beta = vec::dot(d,z);

		x = x+(alpha/beta)*d;
		g = g+(alpha/beta)*z;
		s = M*g;

		beta = alpha;
		alpha = vec::dot(g,s);

		d = (alpha/beta)*d-s;
		l2 = vec::norm(g);
		if (l2<tolerance) {converged=true;break;}
	}

	if (!converged)	cerr<<"PCG failed to converge, norm(g) = "<<l2<<endl;
    return converged;
}

/*solves non-linear Poisson equation using Gauss-Seidel*/
bool PotentialSolver::solveGS()
{
    //references to avoid having to write world.phi
	Field &phi = world.phi;
    Field &rho = world.rho;		//rho contains only ion contribution

	//precompute 1/(dx^2)
    double3 dh = world.getDh();
    double idx2 = 1.0/(dh[0]*dh[0]);
    double idy2 = 1.0/(dh[1]*dh[1]);
    double idz2 = 1.0/(dh[2]*dh[2]);

    double L2=0;			//norm
    bool converged= false;

    /*solve potential*/
    for (unsigned it=0;it<max_solver_it;it++)
    {
 		for (int i=0;i<world.ni;i++)
        	for (int j=0;j<world.nj;j++)
            	for (int k=0;k<world.nk;k++)
                {
                    /*skip over solid (fixed) nodes = Dirichlet boundaries*/
                    if (world.object_id[i][j][k]>0) continue;

                    if (i==0)
                    	phi[i][j][k] = phi[i+1][j][k];
                    else if (i==world.ni-1)
                    	phi[i][j][k] = phi[i-1][j][k];
                    else if (j==0)
                    	phi[i][j][k] = phi[i][j+1][k];
                    else if (j==world.nj-1)
                    	phi[i][j][k] = phi[i][j-1][k];
                    else if (k==0)
                    	phi[i][j][k] = phi[i][j][k+1];
                    else if (k==world.nk-1)
                    	phi[i][j][k] = phi[i][j][k-1];
                    else {	//standard internal open node

                    	//evaluate electron density from the Boltzmann relationshp
                    	double ne = n0 * exp((phi[i][j][k]-phi0)/Te0);

                    	double phi_new = ((rho[i][j][k]-Const::QE*ne)/Const::EPS_0 +
                                        idx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
                                        idy2*(phi[i][j-1][k]+phi[i][j+1][k]) +
                                        idz2*(phi[i][j][k-1]+phi[i][j][k+1]))/(2*idx2+2*idy2+2*idz2);

                    	/*SOR*/
                    	phi[i][j][k] = phi[i][j][k] + 1.4*(phi_new-phi[i][j][k]);
                    }
                }

		 /*check for convergence*/
		 if (it%25==0)
		 {
			double sum = 0;
			for (int i=0;i<world.ni;i++)
				for (int j=0;j<world.nj;j++)
					for (int k=0;k<world.nk;k++)
					{
						/*skip over solid (fixed) nodes*/
						if (world.object_id[i][j][k]>0) continue;

						double R = 0;
                    	if (i==0)
                    		R = phi[i][j][k] - phi[i+1][j][k];
                    	else if (i==world.ni-1)
                    		R = phi[i][j][k] - phi[i-1][j][k];
                    	else if (j==0)
                    		R = phi[i][j][k] - phi[i][j+1][k];
                    	else if (j==world.nj-1)
                    		R = phi[i][j][k] - phi[i][j-1][k];
                    	else if (k==0)
                    		R = phi[i][j][k] - phi[i][j][k+1];
                    	else if (k==world.nk-1)
                    		R = phi[i][j][k] - phi[i][j][k-1];
                    	else {
                    			//evaluate electron density from the Boltzmann relationshp
                    		    double ne = n0 * exp((phi[i][j][k]-phi0)/Te0);
                    			R = -phi[i][j][k]*(2*idx2+2*idy2+2*idz2) +
									(rho[i][j][k]-Const::QE*ne)/Const::EPS_0 +
									idx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
									idy2*(phi[i][j-1][k]+phi[i][j+1][k]) +
									idz2*(phi[i][j][k-1]+phi[i][j][k+1]);
						}

						sum += R*R;
					}

			L2 = sqrt(sum/(world.ni*world.nj*world.nk));
			if (L2<tolerance) {converged=true;break;}
		}
    }

    if (!converged) cerr<<"GS failed to converge, L2="<<L2<<endl;
    return converged;
}

/*solves non-linear Poisson equation using Gauss-Seidel*/
bool PotentialSolver::solveGSLinear(Matrix &A, dvector &x, dvector &b)
{
    double L2=0;			//norm
    bool converged= false;

    /*solve potential*/
    for (unsigned it=0;it<max_solver_it;it++)
    {
 		for (int u=0;u<A.nu;u++)
		{
			double S = A.multRow(u,x)-A(u,u)*x[u]; //multiplication of non-diagonal terms
 			double phi_new = (b[u]- S)/A(u,u);

 			/*SOR*/
            x[u] = x[u] + 1.*(phi_new-x[u]);
		}

		 /*check for convergence*/
		 if (it%25==0)
		 {
			 dvector R = A*x-b;
			 L2 = vec::norm(R);
			 if (L2<tolerance) {converged=true;break;}
		}
    }

    if (!converged) cerr<<"GS failed to converge, L2="<<L2<<endl;
    return converged;
}


/*computes electric field = -gradient(phi) using 2nd order differencing*/
void PotentialSolver::computeEF()
{
	//grab references to data
	Field &phi = world.phi;
	Field3 &ef = world.ef;

	double3 dh = world.getDh();
	double dx = dh[0];
	double dy = dh[1];
	double dz = dh[2];

	for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
			for (int k=0;k<world.nk;k++)
			{
				/*x component*/
				if (i==0)
					ef[i][j][k][0] = -(-3*phi[i][j][k]+4*phi[i+1][j][k]-phi[i+2][j][k])/(2*dx);	/*forward*/
				else if (i==world.ni-1)
					ef[i][j][k][0] = -(phi[i-2][j][k]-4*phi[i-1][j][k]+3*phi[i][j][k])/(2*dx);	/*backward*/
				else
					ef[i][j][k][0] = -(phi[i+1][j][k] - phi[i-1][j][k])/(2*dx);	/*central*/

				/*y component*/
				if (j==0)
					ef[i][j][k][1] = -(-3*phi[i][j][k] + 4*phi[i][j+1][k]-phi[i][j+2][k])/(2*dy);
				else if (j==world.nj-1)
					ef[i][j][k][1] = -(phi[i][j-2][k] - 4*phi[i][j-1][k] + 3*phi[i][j][k])/(2*dy);
				else
					ef[i][j][k][1] = -(phi[i][j+1][k] - phi[i][j-1][k])/(2*dy);

				/*z component*/
				if (k==0)
					ef[i][j][k][2] = -(-3*phi[i][j][k] + 4*phi[i][j][k+1]-phi[i][j][k+2])/(2*dz);
				else if (k==world.nk-1)
					ef[i][j][k][2] = -(phi[i][j][k-2] - 4*phi[i][j][k-1]+3*phi[i][j][k])/(2*dz);
				else
					ef[i][j][k][2] = -(phi[i][j][k+1] - phi[i][j][k-1])/(2*dz);
			}
}

////////////////////Output.cpp//////////////////////

/*saves fields in VTK format*/
void Output::fields(World &world, vector<Species> &species)
{
	/*update gas macroscopic properties*/
	for (Species &sp:species)
		sp.computeGasProperties();

	stringstream name;
	name<<"results/fields_"<<setfill('0')<<setw(5)<<world.getTs()<<".vti";

    /*open output file*/
    ofstream out(name.str());
   	if (!out.is_open()) {cerr<<"Could not open "<<name.str()<<endl;return;}

	/*ImageData is vtk format for structured Cartesian meshes*/
	out<<"<VTKFile type=\"ImageData\">\n";
	double3 x0 = world.getX0();
	double3 dh = world.getDh();
	out<<"<ImageData Origin=\""<<x0[0]<<" "<<x0[1]<<" "<<x0[2]<<"\" ";
	out<<"Spacing=\""<<dh[0]<<" "<<dh[1]<<" "<<dh[2]<<"\" ";
	out<<"WholeExtent=\"0 "<<world.ni-1<<" 0 "<<world.nj-1<<" 0 "<<world.nk-1<<"\">\n";
	
	/*output data stored on nodes (point data)*/
	out<<"<PointData>\n";

	/*object id, scalar*/
	out<<"<DataArray Name=\"object_id\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
	out<<world.object_id;
	out<<"</DataArray>\n";

	/*node volumes, scalar*/
	Field f(world.node_vol);
	for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
			for (int k=0;k<world.nk;k++)
				f[i][j][k] = world.node_vol[i][j][k];
	out<<"<DataArray Name=\"NodeVol\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.node_vol;
	out<<"</DataArray>\n";

	/*potential, scalar*/
	out<<"<DataArray Name=\"phi\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.phi;
	out<<"</DataArray>\n";

	/*charge density, scalar*/
	out<<"<DataArray Name=\"rho\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.rho;
	out<<"</DataArray>\n";

	/*species number densities*/
	for (Species &sp:species)
	{
		out<<"<DataArray Name=\"nd."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.den;
		out<<"</DataArray>\n";
	}
	
	/*species averaged number densities*/
	for (Species &sp:species)
	{
		out<<"<DataArray Name=\"nd-ave."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.den_ave;
		out<<"</DataArray>\n";
	}

	/*species stream velocity*/
	for (Species &sp:species)
	{
		out<<"<DataArray Name=\"vel."<<sp.name<<"\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.vel;
		out<<"</DataArray>\n";
	}

	/*species temperature*/
	for (Species &sp:species)
	{
		out<<"<DataArray Name=\"T."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.T;
		out<<"</DataArray>\n";
	}

	/*electric field, 3 component vector*/
	out<<"<DataArray Name=\"ef\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.ef;
	out<<"</DataArray>\n";

	/*close out tags*/
	out<<"</PointData>\n";

	/*cell data*/
	out<<"<CellData>\n";
	/*species temperature*/
	for (Species &sp:species)
	{
		out<<"<DataArray Name=\"mpc."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.mpc;
		out<<"</DataArray>\n";
	}
	out<<"</CellData>\n";

	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
 	out.close();

 	/*clear samples if not at steady state*/
 	if (!world.isSteadyState())
 		for (Species &sp:species) sp.clearSamples();
}

//writes information to the screen
void Output::screenOutput(World &world, vector<Species> &species)
{
	cout<<"ts: "<<world.getTs();
	for (Species &sp:species)
		cout<<setprecision(3)<<"\t "<<sp.name<<":"<<sp.getNp();
	cout<<endl;
}

//file stream handle
namespace Output {
std::ofstream f_diag;
}

/*save runtime diagnostics to a file*/
void Output::diagOutput(World &world, vector<Species> &species)
{
	using namespace Output;	//to get access to f_diag

	//is the file open?
	if (!f_diag.is_open())
	{
		f_diag.open("runtime_diags.csv");
		f_diag<<"ts,time,wall_time";
		for (Species &sp:species)
			f_diag<<",mp_count."<<sp.name<<",real_count."<<sp.name
				  <<",px."<<sp.name<<",py."<<sp.name<<",pz."<<sp.name
			      <<",KE."<<sp.name;
		f_diag<<",PE,E_total"<<endl;
	}

	f_diag<<world.getTs()<<","<<world.getTime();
	f_diag<<","<<world.getWallTime();

	double tot_KE = 0;
	for (Species &sp:species)
	{
		double KE = sp.getKE();	//species kinetic energy
		tot_KE += KE;		//increment total energy
		double3 mom = sp.getMomentum();

		f_diag<<","<<sp.getNp()<<","<<sp.getRealCount()
			  <<","<<mom[0]<<","<<mom[1]<<","<<mom[2]<<","<<KE;
	}

	//write out system potential and total energy
	double PE = world.getPE();
	f_diag<<","<<PE<<","<<(tot_KE+PE);

	f_diag<<"\n";	//use \n to avoid flush to disc
	if (world.getTs()%25==0) f_diag.flush();
}

/*saves particle x-vel data*/
void Output::particles(World &world, vector<Species> &species, int num_parts) {
	/*loop over all species*/

	for (Species &sp:species) {

		//open a phase_sp_it.vtp
		stringstream name;
		name<<"results/parts_"<<sp.name<<"_"<<setfill('0')<<setw(5)<<world.getTs()<<".vtp";

		/*open output file*/
		ofstream out(name.str());
		if (!out.is_open()) {cerr<<"Could not open "<<name.str()<<endl;return;}

		/*build a list of particles to output*/
		double dp = num_parts/(double)sp.getNp();
		double counter = 0;
		vector<Particle*> to_output;
		for (Particle &part : sp.particles)	{
			counter+=dp;
			if (counter>1) //save particle
				{to_output.emplace_back(&part);counter=-1;}
		}

		/*header*/
		out<<"<?xml version=\"1.0\"?>\n";
		out<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
		out<<"<PolyData>\n";
		out<<"<Piece NumberOfPoints=\""<<to_output.size()<<"\" NumberOfVerts=\"0\" NumberOfLines=\"0\" ";
		out<<"NumberOfStrips=\"0\" NumberOfCells=\"0\">\n";

		/*points*/
		out<<"<Points>\n";
		out<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		for (Particle *part: to_output)
			out<<part->pos<<"\n";
		out<<"</DataArray>\n";
		out<<"</Points>\n";

		/*velocities*/
		out<<"<PointData>\n";
		out<<"<DataArray Name=\"vel."<<sp.name<<"\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		for (Particle *part: to_output)
			out<<part->vel<<"\n";
		out<<"</DataArray>\n";
		out<<"</PointData>\n";

		out<<"</Piece>\n";
		out<<"</PolyData>\n";
		out<<"</VTKFile>\n";

		out.close();
	}

};


