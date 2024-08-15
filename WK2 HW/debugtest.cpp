#include <iostream>
#include <fstream>

using namespace std;

int main() {
    
    ofstream out("debug.csv");
    out << "hi," << endl;

    out << ", , " << "there" << endl;

}

