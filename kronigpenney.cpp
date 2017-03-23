#include "kronigpenney.h"


KronigPenney::KronigPenney() {
    m_accuracy = 1e-5;
    m_waveBasisLength = 0;
    m_unperturbedStatesLength = 0;
    m_perturbedStatesLength = 0;
}

void KronigPenney::setUnitCell(std::__cxx11::string CELLFILE) {
    std::string dummystring;
    std::fstream CELL(CELLFILE , std::ios_base::in);

    if(!CELL.good()){
        CELL.close();
        std::string NewFile = CELLFILE;
        setAReal(vec3(1.0,0,0));
        setBReal(vec3(0,1.0,0));
        setCReal(vec3(0,0,1.0));
        setPotential(-1.0);
        std::ofstream FILE(NewFile);
        FILE << "#CELLPARAMETERS [Angstroms] \n";
        FILE << "1\t0\t0\n0\t1\t0\n0\t0\t1\n\nPotential [eV], Location (a,b,c)\n\n-1\t0.0000\t0.0000\t0.0000\n";
        FILE.close();
    }else{
        CELL >> dummystring >> dummystring;
        double x,y,z, potential;
        CELL >> x >> y >> z;
        setAReal(vec3(x,y,z));
        CELL >> x >> y >> z;
        setBReal(vec3(x,y,z));
        CELL >> x >> y >> z;
        setCReal(vec3(x,y,z));
        CELL >> dummystring >> dummystring >> dummystring >> dummystring;
        CELL >> potential;
        setPotential(potential);
        CELL.close();
    }

    double volume = m_aReal.dot(m_bReal.cross(m_cReal));
    setCellVolume(volume);

    vec3 a,b,c;
    a = TWO_PI*m_bReal.cross(m_cReal)/volume;
    b = TWO_PI*m_cReal.cross(m_aReal)/volume;
    c = TWO_PI*m_aReal.cross(m_bReal)/volume;
    setReciprocalSpace(a,b,c);
}

void KronigPenney::setReciprocalSpace(vec3 a, vec3 b, vec3 c) {
    setAResiprocal(a);
    setBResiprocal(b);
    setCResiprocal(c);
}

void KronigPenney::setWaveBasis(std::__cxx11::string BASISFILE, double energyCutOff)
{
    std::fstream BASIS(BASISFILE, std::ios_base::in);
    double x,y,z, cutOffG;
    int limA, limB, limC;

    if(!BASIS.good()){
        BASIS.close();
        std::string NewFile = BASISFILE;
        std::ofstream FILE(NewFile);

        cutOffG = sqrt(2*ELECTRON_MASS*energyCutOff/(HBAR_C*HBAR_C));
        m_beta = -2.0*log(m_accuracy)/(cutOffG*cutOffG);
        limA = int(cutOffG/m_aResiprocal.length());
        limB = int(cutOffG/m_bResiprocal.length());
        limC = int(cutOffG/m_cResiprocal.length());
//        std::cout << limA << "\t" << limB << "\t" << limC<< std::endl;
        vec3 G;

        double energy;
        for(int l = -limC; l<=limC; l++){
            for(int m = -limB; m<=limB; m++){
                for(int n = -limA; n<=limA; n++){
                    G=n*m_aResiprocal+m*m_bResiprocal+l*m_cResiprocal;
                    energy = HBAR_C*HBAR_C*G.lengthSquared()/(2*ELECTRON_MASS);
                    //std::cout << G << "\t" << energy;
                    if(energy < energyCutOff){
                        FILE << G.x() <<"\t" << G.y() << "\t" << G.z() << "\n";
                        m_waveBasis.push_back(G);
                        m_waveBasisLength ++;
                    }
                }
            }
        }

    }else{
        while(BASIS >> x >> y >> z){
            m_waveBasis.push_back(vec3(x,y,z));
            m_waveBasisLength ++;
        }
        BASIS.close();
    }
}

void KronigPenney::setWaveStates(std::__cxx11::string WAVEFILE, vec3 kPoint) {
    std::fstream WAVES(WAVEFILE, std::ios_base::in);

    double energy;
    std::vector<double> weights;

    if(!WAVES.good()){ //If kPoint is not contained in file
        WAVES.close();
        m_unperturbedStates.clear();
        m_unperturbedStatesLength = 0;

        std::string NewFile = WAVEFILE;
        std::ofstream FILE(NewFile);
        vec3 totalVector;
        for(int i = 0; i<m_waveBasisLength; i++) {
            totalVector = kPoint+m_waveBasis[i];
            energy = (HBAR_C*HBAR_C*totalVector.lengthSquared()/(2*ELECTRON_MASS));

            FILE << kPoint.x() << "\t" << kPoint.y() << "\t" << kPoint.z() ;
            FILE << "\t" << energy;

            weights = std::vector<double>(m_waveBasisLength, 0);
            weights[i] = 1;

            for(int j=0; j<m_waveBasisLength; j++) {
                FILE << "\t" << weights[j] ;
            }
            FILE << "\n";

            m_unperturbedStates.push_back(waveState(kPoint,m_waveBasis[i],energy,weights));
            m_unperturbedStatesLength ++;
        }
        FILE.close();
    }else{
        m_unperturbedStates.clear();
        m_unperturbedStatesLength = 0;
        double x,y,z;
        double weight;
        vec3 kPointFromFile, effectiveG;
        while(WAVES) {
            WAVES >> x >> y >> z;
            WAVES >> energy;
            kPointFromFile = vec3(x,y,z);

            weights.clear();
            for(int i = 0; i<m_waveBasisLength; i++) {
                WAVES >> weight;
                weights.push_back(weight);
            }
            if( (kPointFromFile-kPoint).length() == 0) {
                for(int j = 0; j<m_waveBasisLength; j++) {
                    effectiveG += m_waveBasis[j]*weights[j];
                }
                m_unperturbedStates.push_back(waveState(kPoint,effectiveG,energy,weights));
                m_unperturbedStatesLength ++;
            }
        }
        m_unperturbedStates.pop_back();
        m_unperturbedStatesLength --;

        WAVES.close();

//        if(m_unperturbedStates == 0) {
//            //DO SOME FANCY SHIT
//        }
    }

}

double KronigPenney::greens(double energy) {
    std::complex<double> sum, epsilon;
    sum = 0;
    epsilon = 1i*m_accuracy;
    for (int i = 0; i<m_unperturbedStatesLength; i++) {
        //std::cout << exp(-m_beta*m_unperturbedStates[i].getG().lengthSquared()*log(m_accuracy)/2.0) << std::endl;
        sum += exp(-m_beta*m_unperturbedStates[i].getG().lengthSquared()/2.0)/(energy+epsilon-m_unperturbedStates[i].energy());
    }
    return sum.real();
}

void KronigPenney::eigenValues(vec3 kPoint, double energyMin, double energyMax) {
    double energy, energyStep, green, greenOld, greenCriteria;
    greenCriteria = m_cellVolume/m_potential;
    greenOld = -m_accuracy;
    green = 0;
    energy = energyMin;
    energyStep = 0.005;

    setWaveStates("WAVEFILE_OLD",kPoint);
    std::ofstream WAVES("WAVEFILE",std::ios::app);

    vec3 effectiveG;
    double eigenEnergy, weight, normalizedWeight, basisWeight, normalization;
    std::vector<double> weights, normalizedWeights, basisWeights;
    while (energy <= energyMax) {
        green = greens(energy);
        if (green <= greenCriteria && greenOld >= greenCriteria) {
            std::cout << "FUCCING WANKER! I got some energy!" << std::endl;
            eigenEnergy = energy-(green-greenCriteria)/(green-greenOld)*energyStep;
            WAVES << kPoint.x() << "\t" << kPoint.y() << "\t" << kPoint.z() << "\t" << eigenEnergy;

//            //FIND WAVEFUNCTION
            for(int i = 0; i<m_unperturbedStatesLength; i++) {
                weight = exp(-m_beta*m_unperturbedStates[i].getG().lengthSquared()/2.0)/(eigenEnergy-m_unperturbedStates[i].energy());
                weights.push_back(weight);
                normalization += weight*weight;
            }

//            //NORMALIZE
            for(int i = 0; i<m_unperturbedStatesLength; i++) {
                normalizedWeight = weights[i]/sqrt(normalization);
                normalizedWeights.push_back(normalizedWeight);
            }

            for(int i = 0; i<m_waveBasisLength; i++) {
                basisWeight = 0;
                for(int j = 0; j<m_unperturbedStatesLength; j++) {
                    basisWeight += normalizedWeights[j]*m_unperturbedStates[j].weight(i);
                }
                basisWeights.push_back(basisWeight);
                WAVES << "\t" << basisWeight;
            }

            for(int j = 0; j<m_waveBasisLength; j++) {
                effectiveG += m_waveBasis[j]*basisWeights[j];
            }
            m_perturbedStates.push_back(waveState(kPoint,effectiveG,eigenEnergy,basisWeights));

            WAVES << "\n";

//            //CALCULATE ENERGY
//            double potSum = 0;
//            for(int h = -m_planeWavesRange+1; h<m_planeWavesRange; h++){
//                potSum += m_waveFunction[h+(m_planeWavesRange-1)];
//            }

//            double WaveEnergy, eDiff;
//            WaveEnergy = potSum*potSum*m_potential/m_a.length();
//            for(int h = -m_planeWavesRange+1; h<m_planeWavesRange; h++){
//                tempE = (m_initialEnergies[h+(m_planeWavesRange-1)])*m_waveFunction[h+(m_planeWavesRange-1)]*m_waveFunction[h+(m_planeWavesRange-1)];
//                WaveEnergy += tempE;
//            }
//            std::cout << std::endl;
//            eDiff = energy-(green-potentialInverse)/(green-greenOld)*energyStep-WaveEnergy;
//            std::cout  << WaveEnergy << "\t" << energy-(green-potentialInverse)/(green-greenOld)*energyStep << "\t" << eDiff << std::endl;
//            std::cout << std::endl;
        }
        greenOld = green;
        energy += energyStep;
    }

    WAVES.close();
}

vec3 KronigPenney::aReal() const {
    return m_aReal;
}

void KronigPenney::setAReal(const vec3 &aReal) {
    m_aReal = aReal;
}

vec3 KronigPenney::bReal() const {
    return m_bReal;
}

void KronigPenney::setBReal(const vec3 &bReal) {
    m_bReal = bReal;
}

vec3 KronigPenney::cReal() const {
    return m_cReal;
}

void KronigPenney::setCReal(const vec3 &cReal) {
    m_cReal = cReal;
}

double KronigPenney::cellVolume() const {
    return m_cellVolume;
}

void KronigPenney::setCellVolume(double cellVolume) {
    m_cellVolume = cellVolume;
}

vec3 KronigPenney::aResiprocal() const {
    return m_aResiprocal;
}

void KronigPenney::setAResiprocal(const vec3 &aResiprocal) {
    m_aResiprocal = aResiprocal;
}

vec3 KronigPenney::bResiprocal() const {
    return m_bResiprocal;
}

void KronigPenney::setBResiprocal(const vec3 &bResiprocal) {
    m_bResiprocal = bResiprocal;
}

vec3 KronigPenney::cResiprocal() const {
    return m_cResiprocal;
}

void KronigPenney::setCResiprocal(const vec3 &cResiprocal) {
    m_cResiprocal = cResiprocal;
}

double KronigPenney::potential() const
{
    return m_potential;
}

void KronigPenney::setPotential(double potential)
{
    m_potential = potential;
}



