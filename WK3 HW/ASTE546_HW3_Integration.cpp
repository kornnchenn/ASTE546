#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <cmath>
#include <iterator>
using namespace std;


int main() {
    double g = -9.81; //m/s/s
    const double y0 = 1; //meters
    double dt = 0.01; //seconds
    double v, y_analytical, t, v_old, y_old_b, y_old_f, y_backwardE, y_forwardE, v_avg, y_leapfrog, y_old_l, v_old_lf, y_old_lf, v_lf, y_lf;
    double v0 =0;
    double t0 = 0;
    int iterations = 50;

    ofstream out("HW3_Iterations.csv");
    out << "t, v, y backward, y forward, y leapfrog, y leapfrog early, y analytical" << endl;
    
    v_old = v0;
    y_old_b = y0;
    y_old_f = y0;
    y_old_l = y0;
    //initialize a half step back for leapfrog method
    v_old_lf = v0 - (g * (dt/2));
    y_old_lf = y0;

    for (int i = 0; i < iterations; i++) {
        t = static_cast<double>(i) * dt;

        //calculate velocity
        v = v_old + (g * dt);

        //b
        //calulate y with backward Euler 
        y_backwardE = y_old_b + (v * dt);

        //c
        //calculate y with forward Euler
        y_forwardE = y_old_f + (v_old*dt);

        //e
        //calculate leap frog method (averaging every step)
        v_avg = (v + v_old)/2;
        y_leapfrog = y_old_l + (v_avg * dt);

        //f
        //calculate leapfrog with backward Euler with rewinding velocity 
        v_lf = v_old_lf + (g * (dt));
        y_lf = y_old_lf + (v_lf * dt);

        //calculate analytical solution
        y_analytical = 0.5 * (g * t * t) + y0;

        //export values
        out << t << "," << v << "," << y_backwardE << "," << y_forwardE << "," << y_leapfrog << "," << y_lf << "," << y_analytical << endl;

        //set  old to current
        v_old = v;
        y_old_b = y_backwardE;
        y_old_f = y_forwardE;
        y_old_l = y_leapfrog;
        v_old_lf = v_lf;
        y_old_lf = y_lf;
    }


 return 0;
}

