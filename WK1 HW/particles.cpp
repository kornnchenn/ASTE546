#include <math.h>
#include <fstream>

using namespace std;  //to avoid having to write std::

// this is where the code begins
int main() {
	// we need to store x,y,z for particle 1 and particle 2
	// also need u,v,w for particle 1 and particle 2

	// the common types of variables include
	// int for whole numbers
	// float and double for single and double precision floating point values
	// bool for booleans

	// declare double precision variables to hold particle pos and vel
	double x1= 0.00, y1=0.2, z1=0.0;           // contain garbage data
	double x2 = 0.0, y2 = 0.0, z2 = 0.0;
	double u1=20, v1=0, w1=0;
	double u2=0, v2=0, w2=1;
	// always remember to initialize before using

	double m1 = 1;
	double m2 = 100;

	double q1 = -1e-4;  // some arbitrary values
	double q2 = 1e-4;

	double dt = 0.0005;

	constexpr double PI = 3.1415;
	constexpr double EPS_0 = 8.8542e-12;

	ofstream out("output1.csv");
	out<<"ts,x1,y1,z1,x2,y2,z2,u1,v1,w1,u2,v2,w2"<<endl;

	for (int ts = 0; ts<2000; ts++) {

		double r = sqrt((x2-x1)*(x2-x1) +
					    (y2-y1)*(y2-y1) +
		                (z2-z1)*(z2-z1));

		// compute Coulomb force
		double Fx = 1/(4*PI*EPS_0)*q1*q2/(r*r*r)*(x1-x2);
		double Fy = 1/(4*PI*EPS_0)*q1*q2/(r*r*r)*(y1-y2);
		double Fz = 1/(4*PI*EPS_0)*q1*q2/(r*r*r)*(z1-z2);

		if (r<0.001) {
			Fx = Fy = Fz = 0;
		}

		// update particle 1 velocity
		u1 += (Fx/m1)*dt;
		v1 += (Fy/m1)*dt;
		w1 += (Fz/m1)*dt;

		// update particle 2 velocity
		u2 -= (Fx/m2)*dt;
		v2 -= (Fy/m2)*dt;
		w2 -= (Fz/m2)*dt;

		// update particle 1 position
		x1 += u1*dt;
		y1 += v1*dt;
		z1 += w1*dt;

		// update particle 2 position
		x2 = x2 + u2*dt;
		y2 = y2 + v2*dt;
		z2 = z2 + w2*dt;

		out<<ts<<","
			<<x1<<","<<y1<<","<<z1<<","
		    <<x2<<","<<y2<<","<<z2<<","
			<<u1<<","<<v1<<","<<w1<<","
			<<u2<<","<<v2<<","<<w2<<endl;
	}


	return 0;
}



