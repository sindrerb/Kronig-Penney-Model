#ifndef WAVESTATE_H
#define WAVESTATE_H
#include <vector>
#include <iostream>
#include <vec3.h>

class waveState
{
public:
    waveState();
    waveState(int basisLenght);
    waveState(vec3 kPoint, vec3 effectiveG, double energy, std::vector<double> weights);

    std::vector<double> weights() const;
    void setWeights(const std::vector<double> &weights);

    double weight(int basisNumber);
    void setWeight(int basisNumber, const double basisWeight);

    int waveBasisLength() const;
    void setWaveBasisLength(int waveBasisLength);

    double energy() const;
    void setEnergy(double energy);

    vec3 getK() const;
    void setK(const vec3 &value);

    vec3 getG() const;
    void setG(const vec3 &G);

private:
    //WaveBasis
    int m_waveBasisLength;
    std::vector<double> m_weights;

    //EigenEnergy
    double m_energy;
    vec3 m_k;
    vec3 m_G;
};

#endif // WAVESTATE_H
