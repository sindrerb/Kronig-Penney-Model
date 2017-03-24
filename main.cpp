#include <iostream>
#include <kronigpenney.h>
#include <wavestate.h>
#include <vec3.h>
#include <math.h>
#include <fstream>
using namespace std;

int main(int argc, char *argv[])
{
    double energyCutOff = 100; //eV
    string cellFile = "CELLFILE";
    string basisFile = "BASISFILE";
    string waveFile = "WAVEFILE";
    string oldFile = waveFile+"_OLD";
    string resultFile = "RESULTFILE";

    KronigPenney KP;

    int result;

    KP.initializeCELL(cellFile);

    bool basis = KP.readBASISFILE(basisFile);
    if(!basis){
        KP.setWaveBasis(energyCutOff);
        KP.writeBASISFILE(basisFile);
    }

    vec3 k;

    //loop over desired k-points
    for(double h = 0; h<=0.5; h+=0.01){
        k = h*KP.aResiprocal();

        bool states = KP.readWAVEFILE(waveFile,k);
        if(!states){
            KP.setWaveStates(k);
            KP.writeWAVEFILE(waveFile,k);
        }
        KP.calculateEigenValues(k,-10.0,300.0);
        KP.writeRESULTFILE(resultFile);
    }

    return 0;
}
