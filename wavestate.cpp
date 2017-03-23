#include "wavestate.h"

waveState::waveState() {

}

waveState::waveState(int basisLenght) {
    m_waveBasisLength = basisLenght;
    m_weights = std::vector<double>(m_waveBasisLength, 0);
}

waveState::waveState(vec3 kPoint, vec3 effectiveG, double energy, std::vector<double> weights)
{
        m_weights = weights;
        m_k = kPoint;
        m_G = effectiveG;
        m_energy = energy;
}

std::vector<double> waveState::weights() const
{
    return m_weights;
}

void waveState::setWeights(const std::vector<double> &weights)
{
    m_weights = weights;
}

double waveState::weight(int basisNumber) {
    return m_weights[basisNumber];
}

void waveState::setWeight(int basisNumber, const double basisWeight) {
    m_weights[basisNumber] = basisWeight;
}

int waveState::waveBasisLength() const
{
    return m_waveBasisLength;
}

void waveState::setWaveBasisLength(int waveBasisLength)
{
    m_waveBasisLength = waveBasisLength;
}

double waveState::energy() const
{
    return m_energy;
}

void waveState::setEnergy(double energy)
{
    m_energy = energy;
}

vec3 waveState::getK() const
{
    return m_k;
}

void waveState::setK(const vec3 &value)
{
    m_k = value;
}

vec3 waveState::getG() const
{
    return m_G;
}

void waveState::setG(const vec3 &G)
{
    m_G = G;
}
