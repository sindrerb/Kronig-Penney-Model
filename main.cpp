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

    int result;
    KP.setUnitCell(cellFile);
    KP.setWaveBasis(basisFile,20);

    vec3 k;
    k = vec3(0,0,0);
    KP.setWaveStates(waveFile,k);
    k = vec3(0.1,0,0);
    KP.setWaveStates(waveFile,k);

    //Renaming original files
    result = rename("WAVEFILE", "WAVEFILE_OLD");
//    if(result == 0) {
//        cout << "Success!" << endl;
//    }else{
//        cout << "Fail! you fucker!" << endl;
//    }

//    KP.eigenValues(k,-10.0,100.0);

    return 0;
}
