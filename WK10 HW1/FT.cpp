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
        phi_x[i] = (1/ (c * c)) * -sin(c * i);
    }
    
    return phi_x;
}

// Discrete Fourier Transform (DFT)
vector<complex<double>> DFT(vector<complex<double>> function, double L, int N) {
    vector<complex<double>> fs_phi_x(N);
    complex<double> i(0, 1);  // imaginary unit i = sqrt(-1)
    
    for(int k = 0; k<N; k++){
        complex<double> sum = 0.0;
        for(int n = 0; n<N; n++){
            double angle = 2 * M_PI * k * n / N;
            //sum += function[n] * exp(-i * angle);
            fs_phi_x[k] += function[n] * exp(-i * angle);
        }
        //fs_phi_x[k] += sum;
    }

    //incorporate wave number, -k^2 term
    /*
    for (int k = 1; k<N; k++) {
        double k_val = 2 * M_PI * k / L; //wave number
        fs_phi_x[k] = fs_phi_x[k] / (-k_val * k_val);
    }
    fs_phi_x[0] = 0; */
    

    return fs_phi_x;
}

vector<complex<double>> DFT_derivative(vector<complex<double>>& F, double L, int N) {
    vector<complex<double>> DFDT(N);
    vector<double> k(N);
    double dk = 2 * M_PI / L;  
    
    for (int i = 0; i < (N / 2); i++) {
        k[i] = i * dk;
    }

    k[(N / 2)] = 0;

    for (int i = (N / 2) + 1; i < N; i++) {
        k[i] = (i - N) * dk;
    }

    for (int i = 0; i < N; i++) {
        DFDT[i] = F[i] * complex<double>(0, 1) * k[i]; // Apply i*k
    }

    return DFDT;

}


// Inverse Discrete Fourier Transform (IDFT)
vector<complex<double>> I_DFT(const vector<complex<double>>& F) {
    int N = F.size();
    vector<complex<double>> ft_phi_x(N);
    complex<double> i(0, 1);  // imaginary unit i = sqrt(-1)

    for(int n = 0; n < N; n++) {
        complex<double> sum = 0.0;
        for(int k = 0; k < N; k++) {
            double c = 2 * M_PI * k * n / N;
            sum += F[k] * exp(i * c);
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
    vector<complex<double>> fs_phi_x(N);
    vector<complex<double>> DTFT(N);
    vector<complex<double>> ft_phi_x(N);
    
     

    for (int i = 0; i < N; i++) {
        rho[i] = eps0 * sin(c * i);
        function[i] = -rho[i] / eps0;
    }

    //compute the Fourier Transform of the function
    fs_phi_x = DFT(function, L, N);

    //compute the derivative in the transformed domain
    DTFT = DFT_derivative(fs_phi_x, L, N);

    //compute the inverse Fourier Transform of the derivative
    ft_phi_x = I_DFT(DTFT);


    ofstream out("FT.csv");
    out << "x, analytical, FT_real" << endl;
    // Write each vector element to the file in a new line
    for (size_t i = 0; i < N; i++) {
        out << i << "," << analysticalSolution[i] << "," << real(ft_phi_x[i]) << endl;
     }
    out.close();

    
}