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
    waveState(vec3 kPoint, double energy, std::vector<double> weights);

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

private:
    //WaveBasis
    int m_waveBasisLength;
    std::vector<double> m_weights;

    //EigenEnergy
    double m_energy;
    vec3 m_k;
};

#endif // WAVESTATE_H
