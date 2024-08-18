#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;

vector<double> getAnalyticalSolution(double L, int mode, int N){

    vector<double> phi_x(N);
    double c = 2 * M_PI * mode / L;

    for(int i = 0; i<N; i++){
        phi_x[i] = (1/ (c * c)) * sin(c * i);
    }
    
    return phi_x;
}

// Discrete Fourier Transform (DFT)
vector<complex<double>> DFT(vector<complex<double>> function, double L, int mode, int N) {
    double c = 2 * M_PI * mode / L;
    vector<complex<double>> fs_phi_x(N);
    complex<double> i(0, 1);  // imaginary unit i = sqrt(-1)
    
    for(int k = 0; k<N; k++){
        complex<double> sum = 0.0;
        for(int n = 0; n<N; n++){
            double angle = 2 * M_PI * k * n / N;
            sum += function[n] * exp(-i * angle);
        }
        fs_phi_x[k] += sum;
    }

    return fs_phi_x;
}

// Inverse Discrete Fourier Transform (IDFT)
vector<complex<double>> I_DFT(const vector<complex<double>>& F) {
    int N = F.size();
    vector<complex<double>> ft_phi_x(N);
    complex<double> i(0, 1);  // imaginary unit i = sqrt(-1)

    for(int n = 0; n < N; n++) {
        complex<double> sum = 0.0;
        for(int k = 0; k < N; k++) {
            double angle = 2 * M_PI * k * n / N;
            sum += F[k] * exp(i * angle);
        }
        ft_phi_x[n] = sum / static_cast<double>(N);
    }

    return ft_phi_x;
}


int main() {
    double L = 2 * M_PI;
    int mode = 1;
    int N = 64;
    const double eps0 = 8.854187817e-12;
    vector<double> analysticalSolution = getAnalyticalSolution(L, mode, N);


    double c = 2 * M_PI * mode / L;
    vector<complex<double>> rho(N); 
    vector<complex<double>> function(N);
    vector<complex<double>> ft_phi_x(N);
     

    for (int i = 0; i < N; i++) {
        rho[i] = eps0 * sin(c * i);
        function[i] = -rho[i] / eps0;
    }

    vector<complex<double>> fs_phi_x = DFT(function, L, mode, N);
    ft_phi_x = I_DFT(fs_phi_x);

    ofstream out("FT.csv");
    out << "x, analytical, FT_real, FT_imaginary" << endl;
    // Write each vector element to the file in a new line
    for (size_t i = 0; i < N; i++) {
        out << i << "," << analysticalSolution[i] << "," << real(ft_phi_x[i]) << "," << imag(ft_phi_x[i]) << endl;
     }
    out.close();

    
}