#ifndef KRONIGPENNEY_H
#define KRONIGPENNEY_H

#include <string>
#include <vec3.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <wavestate.h>
#include <vector>
#include <complex>

class KronigPenney
{

public:
    KronigPenney();

    //Setting the environment
    void setUnitCell(vec3 a, vec3 b,vec3 c);
    void setReciprocalSpace(vec3 a, vec3 b,vec3 c);
    void readCELLFILE(std::string CELLFILE);
    void writeCELLFILE(std::string CELLFILE);

    //Setting relevant wave basis for the
    void setWaveBasis(double energyCutOff);
    void readBASISFILE(std::string BASISFILE);
    void writeBASISFILE(std::string BASISFILE);

    //Set predefined eigenstates or generate initial states from the basis
    void setWaveStates(vec3 kPoint);
    void readWAVEFILE(std::string WAVEFILE, vec3 kPoint);
    void writeWAVEFILE(std::string WAVEFILE, vec3 kPoint);

    //Calculate new eigenenergies
    double greens(double energy);
    void calculateEigenValues(vec3 kPoint, double energyMin, double energyMax);

    void
    

    vec3 aReal() const;
    void setAReal(const vec3 &aReal);
    
    vec3 bReal() const;
    void setBReal(const vec3 &bReal);
    
    vec3 cReal() const;
    void setCReal(const vec3 &cReal);
    
    double cellVolume() const;
    void setCellVolume(double cellVolume);
    
    vec3 aResiprocal() const;
    void setAResiprocal(const vec3 &aResiprocal);
    
    vec3 bResiprocal() const;
    void setBResiprocal(const vec3 &bResiprocal);
    
    vec3 cResiporcal() const;
    void setCResiporcal(const vec3 &cResiporcal);
    
    vec3 cResiprocal() const;
    void setCResiprocal(const vec3 &cResiprocal);

    double potential() const;
    void setPotential(double potential);

private:
    //Caclulation parameters
    double m_accuracy;
    double m_beta;

    //Real cell parameters
    vec3 m_aReal,m_bReal,m_cReal;
    double m_cellVolume;

    double m_potential;

    //Reciprocal cell parameters
    vec3 m_aResiprocal, m_bResiprocal, m_cResiprocal;

    //WaveBasis
    int m_waveBasisLength;
    std::vector<vec3> m_waveBasis;

    //WaveStates
    int m_unperturbedStatesLength;
    std::vector<waveState> m_unperturbedStates;

    int m_perturbedStatesLength;
    std::vector<waveState> m_perturbedStates;


    //CONSTANTS
    double TWO_PI = M_PI*2;
    double HBAR_C = 1973.0 ;//eVÅ
    double ELECTRON_MASS = 0.511E6; //eV
};

#endif // KRONIGPENNEY_H
