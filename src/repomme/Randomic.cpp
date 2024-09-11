#include "Randomic.hpp"

auto semente = std::chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937_64 twisterengine(semente);

void Randomic::defSeed(double seed)
{
	if (seed != 0)
	{
		std::mt19937_64 twisterengine2(seed);
		twisterengine = twisterengine2;
	}
}

double Randomic::NumeroRandomico(double Min, double Max)
{
	std::uniform_real_distribution<double> dist_real_uniform(Min, Max);
	return dist_real_uniform(twisterengine);
}

int Randomic::NumeroRandomicoInteiro(int Min, int Max)
{
	std::uniform_int_distribution<int> dist_int_uniform(Min, Max - 1);
	return dist_int_uniform(twisterengine);
}
