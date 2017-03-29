#include <iostream>
#include <kronigpenney.h>
#include <wavestate.h>
#include <vec3.h>
#include <math.h>
#include <fstream>
using namespace std;

int main(int argc, char *argv[])
{
    double energyCutOff = 10; //eV
    string cellFile = "CELLFILE";
    string potentialFile = "POTENTIALFILE";
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


    string dummystring;
    fstream POTENTIAL(potentialFile, std::ios_base::in);
    POTENTIAL >> dummystring >> dummystring >> dummystring >> dummystring;

    double potential, x,y,z, volume;
    vec3 a,b,c;
    vec3 kPoint;
    //while(POTENTIAL){ //LOOP OVER THE DEFINED POTENTIALS
        POTENTIAL >> potential >> x >> y >> z;
        cout << "Calculating with a potential of " << potential << " eV." << endl;
        a = x*KP.aReal();
        b = y*KP.bReal();
        c = z*KP.cReal();
        volume = a.dot(b.cross(c));

        //loop over desired k-points
        for(double h = 0; h<=0.5; h+=0.02){
            kPoint = h*KP.aResiprocal();

            bool states = KP.readWAVEFILE(waveFile,kPoint);
            if(!states){
                KP.setWaveStates(kPoint);
                KP.writeWAVEFILE(waveFile,kPoint);
            }
            KP.calculateEigenValues(kPoint,-1,15.0, potential, volume);
            KP.writeRESULTFILE(resultFile);
        }
        rename("WAVEFILE","WAVEFILE_OLD");
        rename("RESULTFILE","WAVEFILE");
        rename("WAVEFILE_OLD","WAVEFILE_ORIGINAL");
    //}
    POTENTIAL.close();

    return 0;
}
