#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <iterator>
using namespace std;

//define class
//specify new data object of type Rnd
class Rnd {
    public:
    //overload the () operator to return random value
    double operator() () {
        return dist(gen);
    }
    
    protected:
    //class member variables , not accessible from outside the class
    std::default_random_engine gen{std::random_device()()}; //must use {} here
    std::uniform_real_distribution<double> dist{0.0,1.0};
};

//function declarations
double mySignalFunction(float t, float r, int window);
double* windowFilterPointer(double* arr, int size, int window);
vector<double> windowFilterReference(vector<double>& vec, int window);
double sumInRangeArray(const double* arr, int size, int i, int window);
double sumInRangeVector(vector<double>& vec, int i, int window);

//main
int main() {
    //print random number
    //cout << rnd()<<endl;

    ofstream out("HW2Part2_SignalProcessor_window3.csv");
    out << "t, ft, filtered ft Array, filtered ft Vector" << endl;

    double ft;
    float t0 =1;
    float tf =5;
    float t;
    float dt;
    int steps = 1001;
    double results[steps];
    double* filteredSignalArray; 
    vector<double> filteredSignalVector;
    int window = 3;
    
    //ran into some more integer division! oops
    dt = ((tf-t0)/static_cast<float>(steps));
    
    //create noisy signal function
    for(int i=0; i<steps; i++){
        Rnd random;
        float r = random();
        t = t0 + (i*dt);
        ft = mySignalFunction(t,r, window);
        results[i] = ft;
        //out << t << "," << ft << endl;
    }

    //create filtered signal function
    filteredSignalArray = windowFilterPointer(results, steps, window);
    //this translates the results array into a vector, in the below code results is a pointer to the beginning of results
    //and the + steps tells it how many spaces to move past the initial location 
    vector<double> vec(results, results + steps);
    filteredSignalVector = windowFilterReference(vec, window);
    
    
    //print results
    for ( int j = 0; j <= steps; j++){
        t = t0 + (j*dt);
        out << t << "," << results[j] << ",";
        out << filteredSignalArray[j] << "," << filteredSignalVector[j] << endl;
    }

    return 0;
}



//function/class definitions
double mySignalFunction(float t, float r, int window){
     double ft;
     //I was not able to define PI using the expression given in the HW handout, I used the one built into math.h instead
     ft = exp(-0.5 * t) * (cos(6* M_PI *t) + 0.5*(-1 + 2*r)*cos(10* M_PI * t));

    return ft;
}

double* windowFilterPointer(double *arr, int size, int window)
{
    //this function passes in a pointer to an array, filters the array and outputs the filtered array of the same size
    //this is more light weight and provides more control 
    //however it is less safe (cause I might mess up), and it has a fixed size and is thus less flexible.
    
    //define and allocate memory for filtered array;
    double* filteredArray = new double[size];

    for (int i = 0; i < size; i++){
        filteredArray[i] = sumInRangeArray(arr, size, i, window) / (window * 2 + 1 );
    }

    return filteredArray;
}

vector<double> windowFilterReference(vector<double>& vec, int window){
    //this function passes in the array as a variable size vector
    //using the vector function provides dynamic sizing and bounds checking
    //which are nice because I'm not that good at C++ and might mess something up :)
    //however it is less efficient :( 
    
    //defines the beginning of a new vector but it's empty, we will iterate adding new points to the end
    //thus it is dynamically defined 
    vector<double> filteredVector;

    for (auto it = vec.begin(); it != vec.end(); it++){

        int i = distance(vec.begin(), it);
        //push_back adds a new value to the end of vector (or adds a new value)
        filteredVector.push_back( sumInRangeVector(vec,i,window) / (window *2 + 1));
    }

    return filteredVector;
}

double sumInRangeArray(const double* arr, int size, int i,int window){
    //i is the index of interest in the input array
    //windowSize is the size of the window of the filter
    //arr is the array
    //size is size of the array 

    double sum = 0;

    int startIndex = max(0,i-window);
    int endIndex = min(size - 1, i + window);

    for(int j = startIndex; j <= endIndex; j++){
        sum += arr[j];
    }

    return sum;
}

double sumInRangeVector(vector<double>& vec, int i, int window){

    double sum = 0;
    int startIndex = max(0,i-window);
    int endIndex = min(static_cast<int>(vec.size()) - 1, i + window);

    for (int j = startIndex; j <= endIndex; j++){
        sum += vec[j];
    }

    return sum;
}