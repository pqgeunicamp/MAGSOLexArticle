#ifndef RANDOMIC_HPP
#define RANDOMIC_HPP

// Bibliotecas padrão C++
#include <random>				// manipulação de números randômicos
#include <chrono>				// manipulação do tempo

class Randomic
{
    public:
		//Variáveis
		
		//Funções
		void defSeed(double);
		double NumeroRandomico(double, double);
		int NumeroRandomicoInteiro(int Min, int Max);
};

#endif // RANDOM_HPP
