#ifndef MG_HPP
#define MG_HPP

// Bibliotecas padrão C++
#include <math.h> 

// Bibliotecas internas
#include "Mistura.hpp"
#include "Constantes.hpp"

namespace CTE = Constantes;

namespace MG
{

    //Funções PO
    bool CH3(int, int, Mistura&);
    bool CH2(int, int, Mistura&);
    bool CH(int, int, Mistura&);
    bool aCH(int, int, Mistura&);
    bool aCfaC(int, int, Mistura&);
    bool aCfnC(int, int, Mistura&);
    bool nCH2(int, int, Mistura&);
    bool nCH(int, int, Mistura&);
    bool aC_CH3(int, int, Mistura&);
    bool aC_CH2(int, int, Mistura&);
    bool aC_CH(int, int, Mistura&);
    bool aC_C(int, int, Mistura&);
    bool aN(int, int, Mistura&);
    bool CH2SH(int, int, Mistura&);
    bool CH2NH2(int, int, Mistura&);
    bool aCqq (int, int, Mistura&);
    bool aCS (int, int, Mistura&);
    bool aCN (int, int, Mistura&);
    bool aCO (int, int, Mistura&);
    bool SH (int, int, Mistura&);
    bool NH2 (int, int, Mistura&);
    bool OH (int, int, Mistura&);
    bool CH2CO (int, int, Mistura&);
    bool CH3CO (int, int, Mistura&);
    bool CHCO (int, int, Mistura&);  
    bool aC_CO (int, int, Mistura&);    
    bool aldeido_CO (int, int, Mistura&);
    bool aldeido_aC_CO (int, int, Mistura&);
    bool C4rio(int, int, Mistura&);
    bool CH2dCH(int, int, Mistura&);
    bool CH2dC(int, int, Mistura&);
    bool CHdCH(int, int, Mistura&);
    bool CHdC(int, int, Mistura&);
    bool CdC(int, int, Mistura&);
    bool CH2dCdCH(int, int, Mistura&);
    bool CH2dCdC(int, int, Mistura&);
    bool CHdCdCH(int, int, Mistura&);
    int AnelPred(int, int, Mistura&);

    //Funções Auxiliares
    std::vector < int > BuscaLigas_e(int, int, Mistura&);
    bool pmul_e(int, int, int, Mistura&);
    bool tal_e(int, int, int, Mistura&);
    bool toel_e(int, int, Mistura&);
    int ContaLigaAnel(int, int, Mistura&);
    int AnelPred(int, int, Mistura&);
    int nrRings_m(int, Mistura&);
    std::vector < int > nrRings_k_m(int, Mistura&);
    std::vector < std::vector < int > > BuscaLigas_Anel(int, Mistura&);
    int nrlt_e(int, int, int, Mistura&);
    bool CH2pertenceKO(int, int, Mistura&);
    bool CHpertenceKO(int, int, Mistura&);
    bool aCpertence_aCCH2aC(int, int, Mistura&);

    std::vector < std::vector < int > > Galdeido_aC_CO(int, Mistura&);
    std::vector < std::vector < int > > GCHCO(int, Mistura&);
    std::vector < std::vector < int > > GCH2COCH2(int, Mistura&);
    std::vector < std::vector < int > > GCH3CO(int, Mistura&);
    std::vector < std::vector < int > > GaC_CO(int, Mistura&);
    std::vector < std::vector < int > > GaCCH2aC(int, Mistura&);

    std::vector < std::vector < int > > GCH2dCH(int, Mistura&);
    std::vector < std::vector < int > > GCH2dC(int, Mistura&);
    std::vector < std::vector < int > > GCHdCH(int, Mistura&);
    std::vector < std::vector < int > > GCHdC(int, Mistura&);
    std::vector < std::vector < int > > GCdC(int, Mistura&);
    std::vector < std::vector < int > > GCH2dCdCH(int, Mistura&);
    std::vector < std::vector < int > > GCH2dCdC(int, Mistura&);
    std::vector < std::vector < int > > GCHdCdCH(int, Mistura&);
};

#endif // MG_HPP