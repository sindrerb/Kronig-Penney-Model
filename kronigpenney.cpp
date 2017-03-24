#include "kronigpenney.h"


KronigPenney::KronigPenney() {
    m_accuracy = 1e-5;
    m_waveBasisLength = 0;
    m_unperturbedStatesLength = 0;
    m_perturbedStatesLength = 0;
}

void KronigPenney::setUnitCell(vec3 a, vec3 b, vec3 c)
{
    setAReal(a);
    setBReal(b);
    setCReal(c);
}

void KronigPenney::setReciprocalSpace(vec3 a, vec3 b, vec3 c) {
    setAResiprocal(a);
    setBResiprocal(b);
    setCResiprocal(c);
}

bool KronigPenney::readCELLFILE(std::__cxx11::string CELLFILE) {
    std::string dummystring;
    std::fstream CELL(CELLFILE , std::ios_base::in);

    if(!CELL.good()){
        std::cout << "Could not find a CELLFILE" << std::endl;
        return false;
    }else{
        CELL >> dummystring >> dummystring;
        double x,y,z, potential;
        vec3 a,b,c;
        CELL >> x >> y >> z;
        a = vec3(x,y,z);
        CELL >> x >> y >> z;
        b = vec3(x,y,z);
        CELL >> x >> y >> z;
        c = vec3(x,y,z);
        CELL >> dummystring >> dummystring >> dummystring >> dummystring;
        CELL >> potential;
        CELL.close();

        setUnitCell(a,b,c);
        setPotential(potential);

        m_cellVolume = m_aReal.dot(m_bReal.cross(m_cReal));
        setCellVolume(m_cellVolume);

        a = TWO_PI*m_bReal.cross(m_cReal)/m_cellVolume;
        b = TWO_PI*m_cReal.cross(m_aReal)/m_cellVolume;
        c = TWO_PI*m_aReal.cross(m_bReal)/m_cellVolume;

        setReciprocalSpace(a,b,c);
        return true;
    }
}

void KronigPenney::writeCELLFILE(std::__cxx11::string CELLFILE) {
    std::string NewFile = CELLFILE;
    std::ofstream FILE(NewFile);
    FILE << "#CELLPARAMETERS [Angstroms] \n";
    FILE << "1\t0\t0\n0\t1\t0\n0\t0\t1\n\nPotential [eV], Location (a,b,c)\n\n-1\t0.0000\t0.0000\t0.0000\n";
    FILE.close();
}

void KronigPenney::initializeCELL(std::__cxx11::string CELLFILE)
{
    std::fstream CELL(CELLFILE , std::ios_base::in);
    if(!CELL.good()){
        CELL.close();
        setUnitCell(vec3(1.0,0,0),vec3(0,1.0,0),vec3(0,0,1.0));
        setPotential(-1.0);
        writeCELLFILE(CELLFILE);
    }else{
        CELL.close();
        readCELLFILE(CELLFILE);
    }
}

void KronigPenney::setWaveBasis(double energyCutOff)
{
    double cutOffG;
    int gLimA, gLimB, gLimC;

    cutOffG = sqrt(2*ELECTRON_MASS_ENERGY*energyCutOff/(HBAR_C*HBAR_C));

    m_beta = -2.0*log(m_accuracy)/(cutOffG*cutOffG);

    gLimA = int(cutOffG/m_aResiprocal.length());
    gLimB = int(cutOffG/m_bResiprocal.length());
    gLimC = int(cutOffG/m_cResiprocal.length());
    vec3 G;

    double energy;
    for(int l = -gLimC; l<=gLimC; l++){
        for(int m = -gLimB; m<=gLimB; m++){
            for(int n = -gLimA; n<=gLimA; n++){
                G=n*m_aResiprocal+m*m_bResiprocal+l*m_cResiprocal;
                energy = HBAR_C*HBAR_C*G.lengthSquared()/(2*ELECTRON_MASS_ENERGY);

                if(energy <= energyCutOff){
                    m_waveBasis.push_back(G);
                    m_waveBasisLength ++;
                }
            }
        }
    }
    std::cout << "Calculated new wave basis from energy cut off" << std::endl;
}

bool KronigPenney::readBASISFILE(std::__cxx11::string BASISFILE) {
    m_waveBasis.clear();
    m_waveBasisLength = 0;
    std::fstream BASIS(BASISFILE, std::ios_base::in);

    if(!BASIS.good()){
        std::cout << "Could not find a BASISFILE" << std::endl;
        return false;
    }else{
        double x,y,z;
        while(BASIS >> x >> y >> z){
            m_waveBasis.push_back(vec3(x,y,z));
            m_waveBasisLength ++;
        }
        BASIS.close();
        return true;
    }
}

void KronigPenney::writeBASISFILE(std::__cxx11::string BASISFILE)
{
    std::fstream BASIS(BASISFILE, std::ios_base::in);
    vec3 G;

    if(!BASIS.good()){
        BASIS.close();
        std::string NewFile = BASISFILE;
        std::ofstream FILE(NewFile);

        for(int i=0; i<m_waveBasisLength; i++){
            G = m_waveBasis[i];
            FILE << G.x() << "\t" << G.y() << "\t" << G.z() << "\n";
        }
        FILE.close();
        std::cout << "Wrote the wave basis to a BASISFILE" << std::endl;

    }else{
        BASIS.close();
        std::cout << "A BASISFILE does already exist." << std::endl;
    }
}

void KronigPenney::setWaveStates(vec3 kPoint) {
    m_unperturbedStates.clear();
    m_unperturbedStatesLength = 0;

    double energy;
    std::vector<double> weights;

    vec3 totalVector;
    for(int i = 0; i<m_waveBasisLength; i++) {
        totalVector = kPoint+m_waveBasis[i];
        energy = (HBAR_C*HBAR_C*totalVector.lengthSquared()/(2*ELECTRON_MASS_ENERGY));

        weights = std::vector<double>(m_waveBasisLength, 0);
        weights[i] = 1;

        m_unperturbedStates.push_back(waveState(kPoint,m_waveBasis[i],energy,weights));
        m_unperturbedStatesLength ++;
    }
    std::cout << "Generated initial states based on wave basis" << std::endl;
}

bool KronigPenney::readWAVEFILE(std::__cxx11::string WAVEFILE, vec3 kPoint) {
    m_unperturbedStates.clear();
    m_unperturbedStatesLength = 0;

    std::fstream WAVES(WAVEFILE, std::ios_base::in);

    double energy;
    std::vector<double> weights;

    if(!WAVES.good()){ //If kPoint is not contained in file
        WAVES.close();
        std::cout << "Could not find a WAVEFILE" << std::endl;
        return false;
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

        //m_unperturbedStates.pop_back();
        if(m_unperturbedStatesLength > 0){
            m_unperturbedStatesLength --;
        }

        if(m_unperturbedStatesLength == 0) {
            std::cout << "The file does not contain the desired k-point" << std::endl;
            return false;
        }else{
            return true;
        }
    }
}

void KronigPenney::writeWAVEFILE(std::__cxx11::string WAVEFILE, vec3 kPoint) {
    std::string NewFile = WAVEFILE;
    std::ofstream FILE(NewFile, std::ios_base::app);

    double energy;
    std::vector<double> weights;

    for(int i = 0; i<m_unperturbedStatesLength; i++) {
        energy = m_unperturbedStates[i].energy();
        weights = m_unperturbedStates[i].weights();

        FILE << kPoint.x() << "\t" << kPoint.y() << "\t" << kPoint.z() ;
        FILE << "\t" << energy;

        for(int j=0; j<m_waveBasisLength; j++) {
            FILE << "\t" << weights[j] ;
        }
        FILE << "\n";
    }
    FILE.close();
    std::cout << "Wrote a WAVEFILE" << std::endl;
}



//void KronigPenney::setWaveBasis(std::__cxx11::string BASISFILE, double energyCutOff)
//{
//    std::fstream BASIS(BASISFILE, std::ios_base::in);
//    double x,y,z, cutOffG;
//    int limA, limB, limC;

//    if(!BASIS.good()){
//        BASIS.close();
//        std::string NewFile = BASISFILE;
//        std::ofstream FILE(NewFile);

//        cutOffG = sqrt(2*ELECTRON_MASS_ENERGY*energyCutOff/(HBAR_C*HBAR_C));
//        m_beta = -2.0*log(m_accuracy)/(cutOffG*cutOffG);
//        limA = int(cutOffG/m_aResiprocal.length());
//        limB = int(cutOffG/m_bResiprocal.length());
//        limC = int(cutOffG/m_cResiprocal.length());
////        std::cout << limA << "\t" << limB << "\t" << limC<< std::endl;
//        vec3 G;

//        double energy;
//        for(int l = -limC; l<=limC; l++){
//            for(int m = -limB; m<=limB; m++){
//                for(int n = -limA; n<=limA; n++){
//                    G=n*m_aResiprocal+m*m_bResiprocal+l*m_cResiprocal;
//                    energy = HBAR_C*HBAR_C*G.lengthSquared()/(2*ELECTRON_MASS_ENERGY);
//                    //std::cout << G << "\t" << energy;
//                    if(energy < energyCutOff){
//                        FILE << G.x() <<"\t" << G.y() << "\t" << G.z() << "\n";
//                        m_waveBasis.push_back(G);
//                        m_waveBasisLength ++;
//                    }
//                }
//            }
//        }

//    }else{
//        while(BASIS >> x >> y >> z){
//            m_waveBasis.push_back(vec3(x,y,z));
//            m_waveBasisLength ++;
//        }
//        BASIS.close();
//    }
//}

//void KronigPenney::setWaveStates(std::__cxx11::string WAVEFILE, vec3 kPoint) {
//    std::fstream WAVES(WAVEFILE, std::ios_base::in);

//    double energy;
//    std::vector<double> weights;

//    if(!WAVES.good()){ //If kPoint is not contained in file
//        WAVES.close();
//        m_unperturbedStates.clear();
//        m_unperturbedStatesLength = 0;

//        std::string NewFile = WAVEFILE;
//        std::ofstream FILE(NewFile);
//        vec3 totalVector;
//        for(int i = 0; i<m_waveBasisLength; i++) {
//            totalVector = kPoint+m_waveBasis[i];
//            energy = (HBAR_C*HBAR_C*totalVector.lengthSquared()/(2*ELECTRON_MASS_ENERGY));

//            FILE << kPoint.x() << "\t" << kPoint.y() << "\t" << kPoint.z() ;
//            FILE << "\t" << energy;

//            weights = std::vector<double>(m_waveBasisLength, 0);
//            weights[i] = 1;

//            for(int j=0; j<m_waveBasisLength; j++) {
//                FILE << "\t" << weights[j] ;
//            }
//            FILE << "\n";

//            m_unperturbedStates.push_back(waveState(kPoint,m_waveBasis[i],energy,weights));
//            m_unperturbedStatesLength ++;
//        }
//        FILE.close();
//    }else{
//        m_unperturbedStates.clear();
//        m_unperturbedStatesLength = 0;
//        double x,y,z;
//        double weight;
//        vec3 kPointFromFile, effectiveG;
//        while(WAVES) {
//            WAVES >> x >> y >> z;
//            WAVES >> energy;
//            kPointFromFile = vec3(x,y,z);

//            weights.clear();
//            for(int i = 0; i<m_waveBasisLength; i++) {
//                WAVES >> weight;
//                weights.push_back(weight);
//            }
//            if( (kPointFromFile-kPoint).length() == 0) {
//                for(int j = 0; j<m_waveBasisLength; j++) {
//                    effectiveG += m_waveBasis[j]*weights[j];
//                }
//                m_unperturbedStates.push_back(waveState(kPoint,effectiveG,energy,weights));
//                m_unperturbedStatesLength ++;
//            }
//        }
//        m_unperturbedStates.pop_back();
//        m_unperturbedStatesLength --;

//        WAVES.close();

////        if(m_unperturbedStates == 0) {
////            //DO SOME FANCY SHIT
////        }
//    }

//}

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

void KronigPenney::calculateEigenValues(vec3 kPoint, double energyMin, double energyMax) {
    m_perturbedStates.clear();
    m_perturbedStatesLength = 0;
    double energy, energyStep, green, greenOld, greenCriteria;
    greenCriteria = m_cellVolume/m_potential;
    greenOld = -m_accuracy;
    green = 0;
    energy = energyMin;
    energyStep = 0.001;

    double eigenEnergy;
    while (energy <= energyMax) {
        green = greens(energy);
        if (green <= greenCriteria && greenOld >= greenCriteria) {
            eigenEnergy = energy-(green-greenCriteria)/(green-greenOld)*energyStep;
            findPerturbedStates(eigenEnergy, kPoint);
        }
        greenOld = green;
        energy += energyStep;
    }
}

void KronigPenney::findPerturbedStates(double eigenEnergy, vec3 kPoint) {
    vec3 effectiveG;
    std::complex<double> weight, epsilon;
    epsilon = 1i*m_accuracy;
    double normalizedWeight, basisWeight, normalization;
    std::vector<double> weights, normalizedWeights, basisWeights;

    for(int i = 0; i<m_unperturbedStatesLength; i++) {
        weight = exp(-m_beta*m_unperturbedStates[i].getG().lengthSquared()/2.0)/(eigenEnergy+epsilon-m_unperturbedStates[i].energy());
        weights.push_back(weight.real());
        normalization += weight.real()*weight.real();
    }
    //NORMALIZE
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
    }

    for(int j = 0; j<m_waveBasisLength; j++) {
        effectiveG += m_waveBasis[j]*basisWeights[j];
    }
    m_perturbedStates.push_back(waveState(kPoint,effectiveG,eigenEnergy,basisWeights));
    m_perturbedStatesLength ++;
}

void KronigPenney::writeRESULTFILE(std::__cxx11::string RESULTFILE){
    std::ofstream FILE(RESULTFILE, std::ios_base::app);

    double energy;
    std::vector<double> weights;
    vec3 kPoint;
    for(int i = 0; i<m_perturbedStatesLength; i++) {
        energy = m_perturbedStates[i].energy();
        weights = m_perturbedStates[i].weights();
        kPoint = m_perturbedStates[i].getK();

        FILE << kPoint.x() << "\t" << kPoint.y() << "\t" << kPoint.z() ;
        FILE << "\t" << energy;

        for(int j=0; j<m_waveBasisLength; j++) {
            FILE << "\t" << weights[j] ;
        }
        FILE << "\n";
    }
    FILE.close();
    std::cout << "Wrote a new WAVEFILE with perturbed states" << std::endl;
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



