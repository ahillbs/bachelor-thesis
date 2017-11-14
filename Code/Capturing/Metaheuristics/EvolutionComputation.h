//
// Created by gharg on 22.10.17.
//

#ifndef PROJECT_EVOLUTIONCOMPUTATION_H
#define PROJECT_EVOLUTIONCOMPUTATION_H

#include "../Header/CapturingSolution.h"
#include "../Header/Capturing2D.h"
#include "SimulatedAnnealing.h"

using Solutions = std::vector<CapturingSolution>;

namespace Metaheuristic {
	struct EvolveSettings {
		ulong maxTime = 900;
		uint startPopulation = 50;
		uint populationSize = 50;
		ulong seed = 23874697;
		ulong maxIteration = ulong() - 1;
		bool useIteratedMutation = false;
		ulong maxNoChangedBase = 100;
		SASettings innerSASetting;
	};

	class EvolutionComputation {
		EvolveSettings settings;
	public:
		EvolutionComputation();
		EvolutionComputation(EvolveSettings settings);
		const EvolveSettings &getSettings() const;

		void setSettings(const EvolveSettings &settings);

		CapturingSolution solve(Capturing2D *problem);

	private:

		void generatePopulation(std::vector<CapturingSolution> *population, Capturing2D &problem);

		void recombination(Solutions &population, Solutions *nextGeneration);

		void mutateIterated(Solutions *population);

		bool select(Solutions *population, Solutions *secondGen);

		CapturingSolution getBestSolution(Solutions &population);

		void mutateLocalRandom(Solutions *population, Capturing2D *problem);
	};
}

#endif //PROJECT_EVOLUTIONCOMPUTATION_H
