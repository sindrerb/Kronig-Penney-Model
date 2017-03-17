#include "kronigpenney.h"


KronigPenney::KronigPenney() {
    m_waveBasisLength = 0;
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
        std::ofstream FILE(NewFile);
        FILE << "#CELLPARAMETERS [Angstroms] \n";
        FILE << "1\t0\t0\n0\t1\t0\n0\t0\t1\n\nPotential [eV], Location (a,b,c)\n\n-1\t0.0000\t0.0000\t0.0000\n";
        FILE.close();
    }else{
        CELL >> dummystring >> dummystring;
        double x,y,z;
        CELL >> x >> y >> z;
        setAReal(vec3(x,y,z));
        CELL >> x >> y >> z;
        setBReal(vec3(x,y,z));
        CELL >> x >> y >> z;
        setCReal(vec3(x,y,z));
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
        limA = int(cutOffG/m_aResiprocal.length());
        limB = int(cutOffG/m_bResiprocal.length());
        limC = int(cutOffG/m_cResiprocal.length());
        std::cout << limA << "\t" << limB << "\t" << limC<< std::endl;
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

//void KronigPenney::setWaveStates(std::__cxx11::string WAVEFILE,vec3 kPoint) {
//    std::fstream WAVES(WAVEFILE, std::ios_base::in);

//    if(!WAVES.good()){
//        std::cout << "yes" << std::endl;
//        WAVES.close();
//        std::vector<double> weights;
//        weights = std::vector<double>(m_waveBasisLength, 0.0);

//    }else{
//        std::cout << "no" << std::endl;
//        WAVES.close();
//    }
//}

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



