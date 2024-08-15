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

	//Part2 of HW1
	//double x1= 0.00, y1=0.02, z1=0.0;           // contain garbage data
	//double x2 = 1.0, y2 = 0.0, z2 = 0.0;
	//double u1=50, v1=0, w1=0;
	//double u2=-50, v2=0, w2=0;
	//double m1 = 1;
	//double m2 = 10;
	//double q1 = 1e-4;  // some arbitrary values
	//double q2 = 1e-4;

//initialize stuff
	double dt = 0.00005;

	//Part 3:
	constexpr int N = 5;
	double x[N] = {0,0.816,-0.816,0,0.2};
	double y[N] = {0,-0.333,-0.333,0.666,0};
	double z[N] = {0.816,-0.333,-0.333,-0.333,0};
	double u[N] = {0,0,0,0,100};
	double v[N] = {0,0,0,0,-15};
	double w[N] = {0,0,0,0,10};
	double m[N] = {1e6,1e6,1e6,1e6,1};
	double q[N] = {5e-3,5e-3,5e-3,5e-3,-1e-4};
	

	constexpr double PI = 3.1415;
	constexpr double EPS_0 = 8.8542e-12;

	ofstream out("part3output1.csv");
	out<<"ts";
    for (int i=0; i<N; i++){
        out<<",x"<<i<<",y"<<i<<",z"<<i;
    }
    out<<"\n";


	for (int ts = 0; ts<2000; ts++) {

		for (int i=0; i<N; i++)
		{
            double Fx = 0;
	        double Fy = 0;
	        double Fz = 0; 

			for(int j=0; j<5; j++)  {

				double r = sqrt((x[j]-x[i])*(x[j]-x[i]) +
							(y[j]-y[i])*(y[j]-y[i]) +
							(z[j]-z[i])*(z[j]-z[i]));

				// compute Coulomb force
				

				if (i!=j && r>0.001) {
					Fx += 1/(4*PI*EPS_0)*q[i]*q[j]/(r*r*r)*(x[i]-x[j]);
					Fy += 1/(4*PI*EPS_0)*q[i]*q[j]/(r*r*r)*(y[i]-y[j]);
					Fz += 1/(4*PI*EPS_0)*q[i]*q[j]/(r*r*r)*(z[i]-z[j]);
				}
				
				

			}
            // update particle 1 velocity
				u[i] += (Fx/m[i])*dt;
				v[i] += (Fy/m[i])*dt;
				w[i] += (Fz/m[i])*dt;

				// update particle 2 velocity
				//u[j] -= (Fx/m[j])*dt;
				//v[j] -= (Fy/m[j])*dt;
				//w[j] -= (Fz/m[j])*dt;

				// update particle 1 position
				x[i] += u[i]*dt;
				y[i] += v[i]*dt;
				z[i] += w[i]*dt;

				// update particle 2 position
				//x[j] = x[j] + u[j]*dt;
				//y[j] = y[j] + v[j]*dt;
				//z[j] = z[j] + w[j]*dt;

		}
        out<<ts<<","
			    <<x[0]<<","<<y[0]<<","<<z[0]<<","
                <<x[1]<<","<<y[1]<<","<<z[1]<<","
                <<x[2]<<","<<y[2]<<","<<z[2]<<","
                <<x[3]<<","<<y[3]<<","<<z[3]<<","
                <<x[4]<<","<<y[4]<<","<<z[4]<<endl;
		        
	}


	return 0;
}



