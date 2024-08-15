
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

int main() {
    //a)
    float f = 10/4;
    cout << f << endl;

    //output = 2

    float g = 10/4.0;
    cout << g <<endl;

    //output = 2.5

    float h = 10./4;
    cout << h << endl;

    //output = 2.5

    float i = (float)10/4;
    cout << i << endl;

    //output = 2.5

    int ni = 21;
    cout <<1/(ni-1)<<endl;
    cout<<1/(ni-1.0)<<endl;

    //output = 0
    //output = 0.05

    cout <<10%4<<endl;

    //output = 2
 return 0;
}

