#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

struct vector3 { 
    double x, y, z;
};

// Constants
const double E_x = 100.0; //V/m
const double E_y = 0.0; //V/m
const double E_z = 0.0; //V/m

const double B_x = 0.00; // T
const double B_y = 0.00; // T
const double B_z = 0.01; // T

// Particle properties
const double charge = -1.602e-19; // C
const double mass = 3.2e-31;  // kg
const double dt = 1e-12;   //s     
const int iter = 6000; 

void exportToCSV(vector<double> x, vector<double> y) {   
    ofstream out("particle_trajectory.csv");
    out << "x,y\n";
    for (int i = 0; i < iter; i++) {
        out << x[i] << "," << y[i] << "\n";
    }

}

void outputVTP(vector<double> x, vector<double> y,vector<double> u, vector<double> v, int num_points) {
    ofstream out("trace.vtp");
    if (!out.is_open()) {
        cout << "Error opening file" << endl;
        return;
    }

    // //header
    // out << "<?xml version=\"1.0\"?>\n";
    // out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    // out << "<PolyData>\n";
    // out << "<Piece NumberOfPoints=\"" << num_points << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"" << 1 << "\" ";
    // out << "NumberOfStrips=\"0\" NumberOfCells=\"0\">\n";

    // //points
    // out << "<Points>\n";
    // out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    // for (int i = 0; i < num_points; i++) {
    //     out << x[i] << " " << y[i] << " 0";}
    // out <<"\n/DataArray>\n";
    // out << "</Points>\n";

    // //lines
    // out << "<Lines>\n";
    // out << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    // for (int i = 0; i < num_points; i++) {
    //     out << i << " ";}
    // out << "\n/DataArray>\n";
    // out << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    // out << num_points << "\n";
    // out << "</DataArray>\n";
    // out << "</Lines>\n";

    // //point data
    // out << "<PointData>\n";
    // out << "<DataArray Name=\"vel\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    // for (int i = 0; i < num_points; i++) {
    //     out << u[i] << " " << v[i] << " 0 ";}   
    // out << "\n/DataArray>\n";
    
    // out << "<DataArray Name=\"ts\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    // for (int i = 0; i < num_points; i++) out<<i<<" ";
    // out << "\n/DataArray>\n";

    // out << "</PointData>\n";
    // out << "</Piece>\n";
    // out << "</PolyData>\n";
    // out << "</VTKFile>\n";

    // Header
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "<PolyData>\n";
    out << "<Piece NumberOfPoints=\"" << num_points << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << 1<<"\" ";
    out << "NumberOfStrip=\"0\" NumberOfCells=\"0\">\n";

    // Points
    out << "<Points>\n";
    out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < num_points; ++i) {
        out << x[i] << " " << y[i] << " 0\n";
    }
    out << "</DataArray>\n";
    out << "</Points>\n";

    // Lines
    out << "<Lines>\n";
    out << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (int i = 0; i < num_points - 1; ++i) {
        out << i << " " << (i + 1) << "\n";
    }
    out << "</DataArray>\n";
    out << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    out << num_points << "\n";
    out << "</DataArray>\n";
    out << "</Lines>\n";

    // Point data
    out << "<PointData>\n";
    out << "<DataArray Name=\"vel\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < num_points; ++i) {
        out << u[i] << " " << v[i] << " 0\n";
    }
    out << "</DataArray>\n";
    
    out << "<DataArray Name=\"ts\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (int i = 0; i < num_points; ++i) {
        out << i << "\n";
    }
    out << "</DataArray>\n";

    out << "</PointData>\n";
    out << "</Piece>\n";
    out << "</PolyData>\n";
    out << "</VTKFile>\n";

}


int main() {

    vector<double> x(iter);
    vector<double> y(iter);
    vector<double> z(iter);
    vector<double> vx(iter);
    vector<double> vy(iter);
    vector<double> vz(iter);

    x[0] = 0.0;
    y[0] = 0.0;
    z[0] = 0.0;

    vx[0] = 0.0;
    vy[0] = 1.0e5;
    vz[0] = 0.0;


    for (int i = 0; i < iter-1; i++) {
        double Fx, Fy, Fz; 
        Fx = charge * (E_x + (vy[i] * B_z) - (vz[i] * B_y));
        Fy = charge * (E_y + (vz[i] * B_x) - (vx[i] * B_z));
        Fz = charge * (E_z + (vx[i] * B_y) - (vy[i] * B_x));
        double ax = Fx / mass;
        double ay = Fy / mass;
        double az = Fz / mass;

        //update velocity
        vx[i+1] = vx[i] + (ax * dt);
        vy[i+1] = vy[i] + (ay * dt);
        vz[i+1] = vz[i] + (az * dt);

        //update position
        x[i+1] = x[i] + (vx[i] * dt);
        y[i+1] = y[i] + (vy[i] * dt);
        z[i+1] = z[i] + (vz[i] * dt);  
    }
    

    exportToCSV(x, y);

    outputVTP(x, y, vx, vy, iter);


    return 0;
}
