#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

int main (){

    ofstream out("HW2Part1Overflow.csv");
	out<<"i, char value, unsigned char, short, unsigned short"<<endl;

    //char
    char c;
    unsigned char uc;
    short s;
    uc =0;
    c = 0;
    s=0;
    unsigned short us;

    for (int i=0; i<100000; i++) {
        c = i;
        uc = i;
        s =i;
        us =i;
        out<< i <<"," << (int)c << "," << (int)uc << "," << s <<"," << us << endl;
    }
 return 0;
}