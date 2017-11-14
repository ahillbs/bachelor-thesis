//
// Created by gharg on 19.09.17.
//

#ifndef PROJECT_SOLVERWRAPPER_H
#define PROJECT_SOLVERWRAPPER_H

#include "../Header/Capturing2D.h"

class Capturing2D_IP_solver {
	bool cplexOutput = false;
	bool preSolve = false;
public:
	bool isPreSolve() const {
		return preSolve;
	}

	void setPreSolve(bool preSolve) {
		Capturing2D_IP_solver::preSolve = preSolve;
	}
	
	bool isCplexOutput() const {
		return cplexOutput;
	}

	void setCplexOutput(bool cplexOutput) {
		Capturing2D_IP_solver::cplexOutput = cplexOutput;
	}

	CapturingSolution solve(Capturing2D* p) {
		p->setCplexOutput(cplexOutput);
		if(preSolve) {
			CapturingSolution preSol = p->solveGreedy();
			std::cout << "Pre solved: " << preSol.getSolutionValue() << " points: ";
			std::cout << std::endl;
			p->setCplexPreSolution(preSol);
		}
		auto result = p->solveKCapturing();
		return result;
	}
};

class Capturing2D_Greedy_solver {
public:
	CapturingSolution solve(Capturing2D* p) {
		auto result = p->solveGreedy();
		return result;
	}
};

class Capturing2D_IteratedSearch_solver {
	bool localSearch = false;
public:
	bool isLocalSearch() const {
		return localSearch;
	}

	void setLocalSearch(bool localSearch) {
		Capturing2D_IteratedSearch_solver::localSearch = localSearch;
	}

	CapturingSolution solve(Capturing2D* p) {

		CapturingSolution result = p->solveGreedy();
		return result.iteratedLocalSearch(localSearch);
	}
};

class Capturing2D_Random_solution {
private:
	ulong seed = 23564563425;
public:
	ulong getSeed() const {
		return seed;
	}

	void setSeed(ulong seed) {
		Capturing2D_Random_solution::seed = seed;
	}

	CapturingSolution solve(Capturing2D* p) {

		ulong size = p->getGraph()->getSize();
		ulong k = p->getK();
		std::set<Vertex> vertices;

		std::mt19937_64 random(seed);
		std::uniform_int_distribution<Vertex> distribution(0,size);
		while(vertices.size() < k) {
			vertices.insert(distribution(random));
		}

		return CapturingSolution(vertices,p);
	}
};
#endif //PROJECT_SOLVERWRAPPER_H
