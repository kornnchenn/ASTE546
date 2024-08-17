#include <iostream>
#include <ofstream>
#include <math.h>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;

vector<double> getAnalyticalSolution(double L, int mode, int N){

    vector<double> phi_x;
    double c = 2 * M_PI * mode / L;

    for(int i = 0; i<N; i++){
        phi_x[i] = (1/ (c * c)) * sin(c * i);
    }
    
    return phi_x;
}

vector<complex<double>> DFT(double L, int mode, int N) {

    vector<double> ft_phi_x;
    double c = 2 * M_PI * mode / L;

    for(int i = 0; i<N; i++){
        ft_phi_x[i] = ;
    }

    return ft_phi_x;
}

vector<double> IDFT(const vector<complex<double>>& F) {
    
}


int main() {
    double L = 2 * M_PI;
    int mode = 1;
    int N = 64;
    const double eps0 = 8.854187817e-12;
    vector<double> analysticalSolution;


    ofstream out("FT.csv");
    out << "ts, analytical, FT" << endl;
    // Write each vector element to the file in a new line
    for (size_t i = 0; i < analysticalSolution.size(); ++i) {
        out << analysticalSolution[i] << endl;
     }

    
}