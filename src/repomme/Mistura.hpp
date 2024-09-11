#ifndef MISTURA_HPP
#define MISTURA_HPP

// Bibliotecas padrão C++
#include <iostream>				// criação de Input/output files
#include <fstream>              // std::ifstream
#include <vector>				// cria e manipula vetores
#include <cmath>				// antigo <math.h> - libera funções matemáticas, ex. pow
#include <sstream>

// Bibliotecas internas
#include "Randomic.hpp"
#include "Constantes.hpp"		//

#include <Eigen/Dense>

namespace CTE = Constantes;

class Mistura
{
    public:
		//variaveis

        //Informações da mistura
        int nMol;
        std::vector <CTE::Molecula> mol; //vetor de moléculas
        Eigen::MatrixXd z;     //fração molar
        Eigen::MatrixXd zwt;   //fração mássica
        Eigen::MatrixXd Mjoback_moleculas;
        Eigen::MatrixXd Matomos_moleculas;
        Eigen::MatrixXd T_Matomos_moleculas;
        Eigen::MatrixXd Mprop_moleculas;
        Eigen::MatrixXd Mmw;
        Eigen::MatrixXd MW_media;
        Eigen::MatrixXd SOLexEigen;
        Eigen::MatrixXd MGgroupsPO;
        Eigen::MatrixXd f;
        //Eigen::MatrixXd kij;
        Randomic Rnd;
        //unsigned int corte;

        //Dimensões 
		bool aviso;

        typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MyMatrix;
        // inline std::vector < std::vector < double > > ParMGPO;
        
        //Funções
        Mistura();
        Mistura(int);
        void CarregaMistura(std::string = "./run", std::string = "Composicao.csv", std::string = "SOLex.csv", std::string = "MAG.csv");
        void LeMAG(std::string = "./run", std::string = "MAG.csv");
        void LeSOLex(std::string = "./run", std::string = "SOLex.csv");
        void LeComposicao(std::string = "./run", std::string = "Composicao.csv");
        void CarregaMAGnodes(int = -1);
        void ContaUltimoC(int = -1);
        void defSeed(int);
        int BuscaLiga(int, int, int);
        void CarregaLigacoesPorAtomo(int = -1, bool = true);

        void ImprimeMAG(int = -1);
        void ImprimeSOLex(int = -1);
        void ImprimeComposicao(int = -1);
        void ImprimeMAGnode(int = -1);
        void ImprimeLigacoesPorAtomo(int = -1);
        void CarregaSOLexEigen();
        
        //auxiliares
        std::vector < std::vector < double > > LeMatriz(std::string, bool = false);
};

#endif // MISTURA_HPP