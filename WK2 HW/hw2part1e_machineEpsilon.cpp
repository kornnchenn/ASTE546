#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

int main(){
    float f = 0;
    double d = 0;

    ofstream out("HW2Part1e_MachineEpsilon.csv");
    out<<"i, f,d"<<endl;    

    for (int i =0; i<20000000; i++){
        f += 1.0f;
        d += 1.0;

        if (i%100 ==0){
            out << i <<","<< f << "," << d <<endl;
        }
    }

return 0;
}