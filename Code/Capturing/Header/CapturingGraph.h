//
// Created by gharg on 10.07.17.
//

#ifndef BA_CAPTURINGGRAPH_H
#define BA_CAPTURINGGRAPH_H

#include <set>
#include "../../cplex/cplex.hpp"
#include "../../Utilities/Graph/undirectedAdjGraph.h"
#include "CapturingSolution.h"


class CapturingGraph {


protected:
	undirectedAdjGraph graph;
	std::vector<std::vector<Vertex>> vertexFamilies;
	std::map<Vertex, std::set<Vertex>> vertexInFamiliesMapping;
	unsigned int maxK;

	CapturingSolution capSolution;
	bool CplexOutput = true;
	double cplexTimeLimit = 900;
public:
	CapturingGraph(unsigned  int k);

	std::map<Vertex, std::set<Vertex>> & getVertexInFamiliesMapping();

	void setCplexPreSolution(CapturingSolution &solution);
	void setK(unsigned int k);
	unsigned int getK();
	void setCplexOutput(bool hasOutput);
	virtual CapturingSolution solveKCapturing();
	undirectedAdjGraph *getGraph();
	void setGraph(undirectedAdjGraph &graph);
	std::vector<std::vector<Vertex>> *getVertexFamilies();
	void setVertexFamilies(std::vector<std::vector<Vertex>> &families);
	double getCplexTimeLimit() const;
	void setCplexTimeLimit(double cplexTimeLimit);
	CapturingSolution computeSolution(std::set<Vertex> &solution);
	double getSolutionValue();
	CapturingSolution getPreSolution();


	std::pair<unsigned int, unsigned int> getFamilyIndexPair(unsigned long familyIndex,
															 const std::set<Vertex> &solution);
	double getWeightInFamily(unsigned long familyIndex, Vertex v1, Vertex v2);
	double getWeightInFamilyWithIndex(unsigned long familyIndex, unsigned long index1, unsigned long index2);


};


#endif //BA_CAPTURINGGRAPH_H
