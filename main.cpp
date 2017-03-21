#include <iostream>
#include <kronigpenney.h>
#include <wavestate.h>
#include <vec3.h>
#include <math.h>
#include <fstream>
using namespace std;

int main(int argc, char *argv[])
{
    string cellFile = "CELLFILE";
    string basisFile = "BASISFILE";
    string waveFile = "WAVEFILE";
    KronigPenney KP;

    KP.setUnitCell(cellFile);
    KP.setWaveBasis(basisFile,200);
    KP.setWaveStates(waveFile,vec3(0,0,0));
//    vec3 a,b,c;
//    a = vec3(1,0,0);
//    b = vec3(cos(gamma),sin(gamma),0);
//    c = vec3(cos(beta),0,sin(beta));

//    cout << a.length() << endl;
//    cout << b.length() << endl;
//    cout << c.length() << endl;
//    cout << (a+b).length() << endl;
//    cout << (a-b).length() << endl;
//    cout << (a+c).length() << endl;
//    cout << (a-c).length() << endl;
//    cout << (a+b+c).length() << endl;
//    cout << (a-b+c).length() << endl;
//    cout << (a+b+c).length() << endl;
//    cout << (a-b-c).length() << endl;

    return 0;
}
