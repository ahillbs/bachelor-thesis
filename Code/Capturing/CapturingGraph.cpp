//
// Created by gharg on 10.07.17.
//

#include "Header/CapturingGraph.h"

using EdgeToCplexPair = std::pair<std::pair<Vertex ,Vertex >,int>;
using VertexPair = std::pair<Vertex,Vertex>;
CapturingGraph::CapturingGraph(unsigned int k) :  maxK(k), capSolution(this){}//,env(), model(env), cplex(model){ }

void CapturingGraph::setK(unsigned int k) {
	this->maxK = k;
}

unsigned int CapturingGraph::getK() {
	return this->maxK;
}

using SegmentEdge = std::pair<ulong,std::pair<Vertex,Vertex>>;
using EdgeToCplex =  std::map<SegmentEdge,int>;
using CplexToEdge = std::map<ulong,VertexPair >;
template <typename T, typename U, typename V>
inline std::pair<T,std::pair<U,V>> make_triplet(T t, U u, V v) {
	std::make_pair(t,std::make_pair(u,v));
};

CapturingSolution CapturingGraph::solveKCapturing() {
	//translate graph to cplex
	//generate variables
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	//std::vector<Edge> edgesX = graph.getEdges();
	//auto edgeSize = edgesX.size();
	//std::cout << "edgesize: " << edgeSize << std::endl;
	IloNumVarArray cplexBoolArrayX(env);
	IloNumVarArray cplexBoolArrayY(env);
	//EdgeToCplex edgeToCplexVar;

	/*for(ulong l = 0; l < graph.getSize(); ++l) {
		for (ulong i = l+1; i < graph.getSize(); ++i) {
			auto edge = graph.getEdge(l,i);
			if(edge.second) {
				auto v1 = edge.first.first;
				auto v2 = edge.first.second;
				std::string name = "Edge" + std::to_string(v1) + "," + std::to_string(v2);
				cplexBoolArrayX.add(IloBoolVar(env, name.c_str()));
				edgeToCplexVar.insert(EdgeToCplexPair(VertexPair(v1, v2), counter));
				edgeToCplexVar.insert(EdgeToCplexPair(VertexPair(v2, v1), counter));
				++counter;
			}
		}
	}*/
	/*
	for (unsigned int l = 0; l < edgesX.size(); ++l) {
		std::string name = "Edge" + std::to_string(edgesX[l].first) + "," + std::to_string(edgesX[l].second);
		cplexBoolArrayX.add(IloBoolVar(env,name.c_str()));
		edgeToCplexVar.insert(EdgeToCplexPair(VertexPair(edgesX[l].first,edgesX[l].second),l));
		edgeToCplexVar.insert(EdgeToCplexPair(VertexPair(edgesX[l].second,edgesX[l].first),l));
	}*/
	/*
    std::cout << "edgeToCplexVar: " << std::endl;
    for(auto stvPair : edgeToCplexVar) {
        std::cout << "Edge (" << stvPair.first.first << "," << stvPair.first.second << "): " << stvPair.second << std::endl;
    }*/

	for (unsigned int l = 0; l < graph.getVertices().size(); ++l) {
		std::string name = "Vertex" + std::to_string(l);
		cplexBoolArrayY.add(IloBoolVar(env,name.c_str()));
	}
	model.add(cplexBoolArrayX);
	model.add(cplexBoolArrayY);

	//Adding Constraints
	CplexToEdge cplexToEdge;
	uint familyConstraints = 0;
	uint families = 0;
	uint edgesUsed = 0;
	ulong counter = 0;
	for (auto vertexFamily : vertexFamilies) {

		edgesUsed += vertexFamily.size()-1;
		cplexToEdge.insert(std::make_pair(counter,std::make_pair(vertexFamily[0],vertexFamily[1])));

		std::string name = "Edge" + std::to_string(families) +
				"," +std::to_string(vertexFamily[0]) +
				"," + std::to_string(vertexFamily[1]);
		cplexBoolArrayX.add(IloBoolVar(env, name.c_str()));
		IloRange  rangeStart(env, 0 ,  cplexBoolArrayY[vertexFamily[0]] - cplexBoolArrayX[counter]);
		model.add(rangeStart);
		++familyConstraints;
		familyConstraints += (vertexFamily.size()-2) * 2;
		//edgeUsedCounter[counter]++;
		++counter;
		for (unsigned int i = 1; i < vertexFamily.size()-1; ++i) {


			cplexToEdge.insert(std::make_pair(counter,std::make_pair(vertexFamily[i],vertexFamily[i+1])));
			std::string name = "Edge" + std::to_string(families) +
							   "," +std::to_string(vertexFamily[i]) +
							   "," + std::to_string(vertexFamily[i+1]);
			cplexBoolArrayX.add(IloBoolVar(env, name.c_str()));

			IloRange  range(env,0, cplexBoolArrayY[ vertexFamily[i]] + cplexBoolArrayX[counter] - cplexBoolArrayX[counter-1] );

			IloRange range2(env,0, cplexBoolArrayY[vertexFamily[i]] + cplexBoolArrayX[counter-1] - cplexBoolArrayX[counter] );

			model.add(range);
			model.add(range2);

			//edgeUsedCounter[counter]++;
			++counter;
		}
		//edgeX1 = edgeToCplexVar[g.getEdge(vertexFamily.second[vertexFamily.second.size()-2],vertexFamily.second[vertexFamily.second.size()-1]).first];
		IloRange  rangeEnd(env, 0 ,cplexBoolArrayY[vertexFamily[vertexFamily.size()-1]] - cplexBoolArrayX[counter-1]);
		model.add(rangeEnd);
		++familyConstraints;
		++families;
	}

	// SUM(Vertices) <= maxK Constraint
	IloNumExpr expr2(env);
	for(unsigned int i = 0; i < cplexBoolArrayY.getSize(); ++i) {
		expr2 += cplexBoolArrayY[i];
	}
	model.add(expr2 <= (int)maxK);
	expr2.end();



	counter = 0;
	// Adding maximize Subject
	IloNumExpr maxExpr(env);
	for(unsigned int i = 0; i < cplexBoolArrayX.getSize(); ++i) {
		Edge edge = graph.getEdge(cplexToEdge[i]).first;
		maxExpr += edge.third * cplexBoolArrayX[i];
	}
	model.add(IloMaximize(env, maxExpr));
	maxExpr.end();
	//*/
	///*
	try {
		if (CplexOutput) {
			std::cout << "Amount nodes: " << cplexBoolArrayY.getSize() << " Amount edges: " << cplexBoolArrayX.getSize() << std::endl;
			//cplex.setOut(std::cout);
			std::cout << "Amount constraints: " << familyConstraints << " amount families: " << families << std::endl;
			std::cout << "Amount edges in families: " << edgesUsed << std::endl;

		} else {
			cplex.setOut(env.getNullStream());
		}

		//std::cout << cplex.getModel() << std::endl;
		//if initial solution is present add to cplex as MIPstart
		if(capSolution.getVertexSolution().size() > 0) {
			//fill initial solution
			auto& initVertexSol = capSolution.getVertexSolution();
			auto& initEdgeSol = capSolution.getEdgeSolution();
			IloNumArray startValY(env);
			double capSolutionVal = 0;
			for (int i = 0; i < cplexBoolArrayY.getSize(); ++i) {
				startValY.add(IloBool(initVertexSol.find(i) != initVertexSol.end()));
			}
			std::cout << "Value computed: " << capSolutionVal << std::endl;
			cplex.addMIPStart(cplexBoolArrayY,startValY,IloCplex::MIPStartAuto,"YBoolStart");
		}
		//cplex.setParam(IloCplex::Param::Emphasis::Numerical,true);
		cplex.setParam(IloCplex::Param::TimeLimit,cplexTimeLimit);
		if (cplex.solve()) {
			IloNumArray vals(env);
			if (CplexOutput) {
				env.out() << "Solution status = " <<
						  cplex.getStatus() << std::endl;
				env.out() << "Solution value = " <<
						  cplex.getObjValue() << std::endl;

			}
			cplex.getValues(vals, cplexBoolArrayY);
			std::set<Vertex> vertexSolution;
			std::vector<Edge> edgeSolution;
			for(unsigned int i = 0; i < vals.getSize() ; ++i) {
				if (vals[i] > 0) {
					vertexSolution.insert((Vertex) i);
				}
			}
			cplex.getValues(vals, cplexBoolArrayX);
			for(unsigned int i = 0; i < vals.getSize() ; ++i) {
				if (vals[i] > 0) {
					Edge e = graph.getEdge(cplexToEdge[i]).first;
					edgeSolution.push_back(e);
				}
			}
			double objVal = cplex.getObjValue();
			double upperBound = cplex.getBestObjValue();

			cplex.clear();
			cplex.end();
			model.end();
			env.end();
			CapturingSolution sol(vertexSolution, edgeSolution, objVal, this);
			sol.setUpperValue(upperBound);
			return sol;

		}
	} catch (IloException iloException) {
		iloException.print(std::cout);
	}
	return CapturingSolution(this);
}


undirectedAdjGraph* CapturingGraph::getGraph() {
	return &this->graph;
}

void CapturingGraph::setGraph(undirectedAdjGraph &graph) {
	this->graph = graph;
}

std::vector<std::vector<Vertex>> *CapturingGraph::getVertexFamilies() {
	return &this->vertexFamilies;
}

void CapturingGraph::setVertexFamilies(std::vector<std::vector<Vertex>> &families) {
	this->vertexFamilies = families;
}


void CapturingGraph::setCplexPreSolution(CapturingSolution &solution) {
	/*this->edgeSolution.clear();
	this->vertexSolution = solution;
	for (unsigned int i = 0; i < vertexFamilies.size(); ++i) {

		Vertex startVertex = Vertex() - 1;
		Vertex endVertex = Vertex() - 1;
		unsigned long j = 0;
		unsigned long k = vertexFamilies[i].size() - 1;
		while(k > j && (endVertex == Vertex() - 1 || startVertex == Vertex() - 1)) {
			if(startVertex == Vertex() - 1 && (solution.find(vertexFamilies[i][j]) != solution.end())) {
				startVertex = vertexFamilies[i][j];
			} else if (startVertex == Vertex() - 1) {
				++j;
			}
			if(endVertex == Vertex() - 1 && (solution.find(vertexFamilies[i][k]) != solution.end())) {
				endVertex = vertexFamilies[i][k];
			} else if (endVertex == Vertex() - 1 ) {
				--k;
			}
		}
		auto indexes = getFamilyIndexPair(i,solution);
		unsigned int j = indexes.first;
		unsigned int k = indexes.second;
		//found a segment track
		if(k > j) {
			for (unsigned long l = j+1; l <= k ; ++l) {
				// add every edge to the total length
				Edge edge = graph.getEdge(vertexFamilies[i][l-1],vertexFamilies[i][l]).first;
				capSolution.getEdgeSolution().insert(edge);
			}
		}
	}*/
	this->capSolution = solution;
}

CapturingSolution CapturingGraph::computeSolution(std::set<Vertex> &solution) {
	double length = 0;
	std::vector<Edge> edges;
	for (unsigned int i = 0; i < vertexFamilies.size(); ++i) {
		/*
		Vertex startVertex = Vertex() - 1;
		Vertex endVertex = Vertex() - 1;
		unsigned long j = 0;
		unsigned long k = vertexFamilies[i].size() - 1;
		while(k > j && (endVertex == Vertex() - 1 || startVertex == Vertex() - 1)) {
			if(startVertex == Vertex() - 1 && (solution.find(vertexFamilies[i][j]) != solution.end())) {
				startVertex = vertexFamilies[i][j];
			} else if (startVertex == Vertex() - 1){
				++j;
			}
			if(endVertex == Vertex() - 1 && (solution.find(vertexFamilies[i][k]) != solution.end())) {
				endVertex = vertexFamilies[i][k];
			} else if (endVertex == Vertex() - 1){
				--k;
			}

		}*/
		std::pair<unsigned int, unsigned int> indexes = getFamilyIndexPair(i,solution);
		//found a segment track
		unsigned int j = indexes.first;
		unsigned int k = indexes.second;
		if(k > j) {
			for (unsigned long l = j+1; l <= k ; ++l) {
				// add every edge to the total length
				Edge edge = graph.getEdge(vertexFamilies[i][l-1],vertexFamilies[i][l]).first;
				edges.push_back(edge);
				length += edge.third;
			}
		}
	}
	return CapturingSolution(solution, edges, length, this);
}

double CapturingGraph::getSolutionValue() {
	return this->capSolution.getSolutionValue();
}

/// Returns the two family indexes, which have the longest weight/path
/// \param familyIndex Index of the Family searching
/// \param solution solution to be used in the family
/// \return pair of indexes in the family
std::pair<unsigned int, unsigned int> CapturingGraph::getFamilyIndexPair(unsigned long familyIndex,
																		 const std::set<Vertex> &solution) {
	Vertex startVertex = Vertex() - 1;
	Vertex endVertex = Vertex() - 1;
	unsigned long j = 0;
	unsigned long k = vertexFamilies[familyIndex].size() - 1;
	while(k > j && (endVertex == Vertex() - 1 || startVertex == Vertex() - 1)) {
		if(startVertex == Vertex() - 1 &&
				(solution.find(vertexFamilies[familyIndex][j]) != solution.end())) {
			startVertex = vertexFamilies[familyIndex][j];
		} else if (startVertex == Vertex() - 1) {
			++j;
		}
		if(endVertex == Vertex() - 1 &&
				(solution.find(vertexFamilies[familyIndex][k]) != solution.end())) {
			endVertex = vertexFamilies[familyIndex][k];
		} else if (endVertex == Vertex() - 1 ) {
			--k;
		}
	}
	return std::pair<unsigned int, unsigned int>(j,k);
}

double CapturingGraph::getWeightInFamily(unsigned long familyIndex, Vertex v1, Vertex v2) {
	std::set<Vertex> pairSolution = {v1,v2};
	auto indexes = getFamilyIndexPair(familyIndex,pairSolution);
	unsigned int j = indexes.first;
	unsigned int k = indexes.second;
	if(k > j) {
		double length = 0;
		for (unsigned long l = j+1; l <= k ; ++l) {
			// add every edge to the total length
			Edge edge = graph.getEdge(vertexFamilies[familyIndex][l-1],vertexFamilies[familyIndex][l]).first;
			length += edge.third;
		}
		return length;
	}
	return 0;
}

double CapturingGraph::getWeightInFamilyWithIndex(unsigned long familyIndex, unsigned long index1, unsigned long index2) {
	if (index2 < index1) {
		return 0;
	}
	double length = 0;
	for (unsigned long l = index1 + 1; l <= index2; ++l) {
		// add every edge to the total length
		Edge edge = graph.getEdge(vertexFamilies[familyIndex][l - 1], vertexFamilies[familyIndex][l]).first;
		length += edge.third;
	}
	return length;

}

void CapturingGraph::setCplexOutput(bool hasOutput) {
	this->CplexOutput = hasOutput;
}

CapturingSolution CapturingGraph::getPreSolution() {
	return this->capSolution;
}

std::map<Vertex, std::set<Vertex>> & CapturingGraph::getVertexInFamiliesMapping() {
	return vertexInFamiliesMapping;
}

double CapturingGraph::getCplexTimeLimit() const {
	return cplexTimeLimit;
}

void CapturingGraph::setCplexTimeLimit(double cplexTimeLimit) {
	CapturingGraph::cplexTimeLimit = cplexTimeLimit;
}




