//
// Created by gharg on 6/12/17.
//

#include <random>
#include <functional>
#include <algorithm>
#include <set>
#include <CGAL/box_intersection_d.h>
#include "Header/Capturing2D.h"
#include "../Visualization/SegmentVisualization.h"
template<typename T>
void removeDuplicates(std::vector<T>* vector) {
	auto it = vector->begin();
	T t = *it;
	while(it != vector->end()) {
		++it;
		T t2 = *it;
		if(t == t2) {
			vector->erase(it);
			it = vector->begin();
		}
	}
}

Capturing2D::Capturing2D(unsigned  int k) : CapturingGraph(k), isConverted(false), segmentGraph() { }
Capturing2D::Capturing2D(unsigned int k, std::vector<Segment<NT>> segments, uint xBound = 1000, uint yBound = 1000)
		: CapturingGraph(k), isConverted(false), segmentGraph(), xBound(xBound), yBound(yBound) {
	for (uint i = 0; i < segments.size(); ++i) {
        this->segments.push_back(segments[i]);
		this->segmentGraph.addVertex();
		vectorFamilies.push_back(std::vector<Vector2D>());
		vertexFamilies.push_back(std::vector<Vertex>());
    }
	convertSegments2Graph();
}

//typedef boost::property<boost::edge_weight_t, double> edgeWeightProperty;
//typedef boost::adjacency_list<boost::listS, boost::vecS,boost::undirectedS,boost::no_property,edgeWeightProperty> Graph;
//typedef Graph::vertex_descriptor Vertex;
//typedef Graph::edge_descriptor Edge;

/// Generates random Segments between  0 and xBound and 0 and yBound
/// \param amount Amount of segments to be created
/// \param xBound Bound in x direction
/// \param yBound Bound in y direction
void Capturing2D::generateRandomObjects2D(int amount, uint xBound, uint yBound, ulong seed) {
    this->xBound = xBound;
    this->yBound = yBound;
	this->isConverted = false;
    std::mt19937_64 random(seed);
    std::mt19937_64 randomY(seed+1);
	segmentGraph = undirectedAdjGraph();
	segments.clear();

    auto dice_randX = std::bind(std::uniform_int_distribution<int>(0,xBound), random);
    auto dice_randY = std::bind(std::uniform_int_distribution<int>(0,yBound), randomY);
    for (int i = 0; i < amount; ++i) {
        int coordinates[4];
        for (int j = 0; j < 4; ++j) {

            if(j % 2 == 0) {
                coordinates[j] = dice_randX();
            } else {
                coordinates[j] = dice_randY();
            }
        }
        segments.push_back(Segment<NT>(coordinates[0],coordinates[1],coordinates[2],coordinates[3]));
		segmentGraph.addVertex();
		vectorFamilies.push_back(std::vector<Vector2D>());
		vertexFamilies.push_back(std::vector<Vertex>());
    }
}



void Capturing2D::convertSegments2Graph() {
	std::cout << "start converting" << std::endl;
    std::vector<Vertex> vecs;
    //std::map<unsigned int,std::vector<Vector2D>> vectorFamilies;
	//auto * vecSet = this->calculateSegmentIntersection();

	Segments cgal_segments;
	convertToCGALSegments(segments, &cgal_segments);
	std::set<Point_2> points;
	std::vector<std::set<Point_2>> familyPoints(segments.size());
	std::map<Point_2, std::set<Vertex>> vectorVertexMap;
	computeIntersectionPoints(&cgal_segments, &points, &familyPoints, &vectorVertexMap);
	ulong edgesInFamilies = 0;
	for(auto family : familyPoints) {
		edgesInFamilies += family.size()-1;
	}
	std::set<Vector2D> vecSet;
	convertCGALToGraph(cgal_segments, &points, &familyPoints, &vectorVertexMap, &vecSet);


	// translate points to vertices and edges in Graph
	this->translateGraph(vecSet);

	/*delete(report->segments1);
	delete(report->points);
	delete(report->vectorVertexMap);
	delete(report->familyPoints);
	delete(report->points);
	delete(report);*/
	//delete(vecSet);
	isConverted = true;
}

CapturingSolution Capturing2D::solveKCapturing() {
	if(!isConverted) {
		this->convertSegments2Graph();
	}
	return CapturingGraph::solveKCapturing();
}

void Capturing2D::printGraph() {
	std::string name = "TestGraph";
	printGraph(name);
}

void Capturing2D::printGraph(const char* cName) {
	std::string name(cName);
	printGraph(name);
}

void Capturing2D::printGraph(std::string docName) {
	this->printGraph(docName,capSolution);
}

void Capturing2D::printGraph(const char *cName, const CapturingSolution &solution) {
	std::string docName(cName);
	printGraph(docName, solution);
}

void Capturing2D::printGraph(std::string docName, const CapturingSolution &solution) {
	if(xBound > 0 && yBound > 0  && solution.getGraph() == this) {
		//std::vector<Edge> edgesX = graph.getEdges();
		docName += ".svg";
		SegmentVisualization document(docName.data(), xBound, yBound);
		document.setSegments(vectorFamilies);
		std::set<std::pair<Vector2D,Vector2D>> takenEdges;
		for (uint k = 0; k < vectorFamilies.size(); ++k) {
			auto indexPair = this->getFamilyIndexPair(k, solution.getVertexSolution());
			if (indexPair.first < indexPair.second) {
				takenEdges.insert(std::pair<Vector2D, Vector2D>(vectorFamilies[k][indexPair.first],
																vectorFamilies[k][indexPair.second]));
			}
		}
		/*
		std::set<Edge> edgeSol = solution.getEdgeSolution();
		for (auto i = edgeSol.begin(); i != edgeSol.end(); ++i) {

			auto vec1 = vectorMapping[i->first];
			auto vec2 = vectorMapping[i->second];

			takenEdges.insert(std::pair<Vector2D,Vector2D>(vec1,vec2));

		}
		 */
		std::vector<Vector2D> vertices;
		auto& vertexSol = solution.getVertexSolution();
		for (auto j = vertexSol.begin(); j != vertexSol.end(); ++j) {
			vertices.push_back(vectorMapping[*j]);
		}
		document.setTakenVertices(vertices);
		document.setTakenEdges(takenEdges);
		document.drawSvg();
		document.saveDocument();
		std::cout << "Document " << docName << " saved" << std::endl;
	}
}

void Capturing2D::generateRandomPoints2D(unsigned long pointAmount, unsigned long segmentAmount, uint xBound,
										 uint yBound, ulong seed) {
	this->xBound = xBound;
	this->yBound = yBound;
	this->isConverted = false;
	std::mt19937_64 random(seed);
	std::mt19937_64 randomY(seed+1);
	std::mt19937_64 randomS(seed+2);
	segmentGraph = undirectedAdjGraph();
	segments.clear();

	auto dice_randX = std::bind(std::uniform_int_distribution<int>(0,xBound), random);
	auto dice_randY = std::bind(std::uniform_int_distribution<int>(0,yBound), randomY);
	std::vector<Vector2D> vectors;
	for (unsigned long i = 0; i < pointAmount; ++i) {
		int coordinates[2];
			coordinates[0] = dice_randX();
			coordinates[1] = dice_randY();
			vectors.push_back(Vector2D(coordinates[0],coordinates[1]));
	}
	unsigned long maxEdges = pointAmount*(pointAmount-1) / 2;
	/*std::vector<Segment> allSegments;
	for (unsigned long j = 0; j < maxEdges; ++j) {
		for (unsigned long i = j+1; i < maxEdges; ++i) {
			allSegments.push_back(Segment(vectors[i], vectors[j]));
		}
	}*/

	if(segmentAmount < maxEdges) {
		auto dice_randS = std::bind(std::uniform_int_distribution<unsigned long>(0, vectors.size() - 1), randomS);

		using longPair = std::pair<unsigned long, unsigned long>;
		std::set<longPair> indexPairs;
		for (unsigned long j = 0; j < segmentAmount; ++j) {
			auto vecIndex1 = dice_randS();
			auto vecIndex2 = dice_randS();
			while(vecIndex1 == vecIndex2) {
				vecIndex2 = dice_randS();
			}
			longPair pair = longPair(vecIndex1,vecIndex2);
			if(indexPairs.find(pair) == indexPairs.end()) {
				longPair pair2 = longPair(vecIndex2,vecIndex1);
				indexPairs.insert(pair);
				indexPairs.insert(pair2);
				segmentGraph.addVertex();
				segments.push_back(Segment<NT>(vectors[vecIndex1], vectors[vecIndex2]));
				vectorFamilies.push_back(std::vector<Vector2D>());
				vertexFamilies.push_back(std::vector<Vertex>());
			} else {
				// could be trouble for segmentAmount near maxEdge
				//--j;
			}
		}/*
		unsigned long maxIndex = maxEdges-1;
		auto dice_randS = std::bind(std::uniform_int_distribution<unsigned long>(0, maxIndex*50), randomS);
		for (unsigned long j = 0; j < segmentAmount; ++j) {
			segmentGraph.addVertex();
			auto segIndex = dice_randS();

			segments.push_back(allSegments[segIndex % maxIndex]);
			std::vector<Segment>::const_iterator iterator = allSegments.begin();
			iterator += segIndex % maxIndex;
			allSegments.erase(iterator);
			--maxIndex;
			vectorFamilies.push_back(std::vector<Vector2D>());
			vertexFamilies.push_back(std::vector<Vertex>());
		}*/
	} else {
		std::cout << "Amount of Segments exceed maximum Edges for " << pointAmount << " points! " << std::endl;
		std::cout << "Complete graph will be build instead" << std::endl;
		for (unsigned long j = 0; j < maxEdges; ++j) {
			for (unsigned long i = j+1; i < maxEdges; ++i) {
				segmentGraph.addVertex();
				segments.push_back(Segment<NT>(vectors[i], vectors[j]));
				vectorFamilies.push_back(std::vector<Vector2D>());
				vertexFamilies.push_back(std::vector<Vertex>());
			}
		}
	}

}

CapturingSolution Capturing2D::solveGreedy() {
	return solveGreedy(ulong() - 1);
}

CapturingSolution Capturing2D::solveGreedy(ulong startSegment = ulong() - 1) {
	std::set<Vertex> greedySolution;
	// first we need to construct a subgraph where all vertices have a degree > 2
	std::map<Vector2D, Vertex> vertexMappingSub;
	std::map<Vertex , Vector2D> vectorMappingSub;
	std::vector<std::vector<Vector2D>> vectorFamiliesSub;
	undirectedAdjGraph subGraph = computeConnectedSubgraph(&vertexMappingSub, &vectorMappingSub, &vectorFamiliesSub);

	std::pair<std::pair<Vector2D, Vector2D>, NT> vectorPair;
	if(startSegment != ulong() - 1) {
		vectorPair = getFamilyMaxLength(startSegment,vectorFamiliesSub);
		//if chosen segment will not have any points left after reduce chose best segment instead
		//just a fast workaround
		if (vectorPair.first.first == Vector2D(0,0) && vectorPair.first.second == Vector2D(0,0) && vectorPair.second == 0) {
			vectorPair = getHeaviestFamily(vectorFamiliesSub);
		}
	} else {
		// 1. search longest Segment
		vectorPair = getHeaviestFamily(vectorFamiliesSub);
	}
	std::cout << "heaviest segment found" << std::endl;
	greedySolution.insert(vertexMappingSub[vectorPair.first.first]);
	greedySolution.insert(vertexMappingSub[vectorPair.first.second]);
	NT currentWeight = vectorPair.second;

	// 2. loop which adds the vertex with most weight gain
	while(greedySolution.size() < maxK) {
		std::cout << "iterating..." << greedySolution.size() << std::endl;
		// sort out possible Vertices
		auto bestVertexPair = getBestVertex(greedySolution, vertexMappingSub, vectorMappingSub, vectorFamiliesSub);

		//if no other connection is found there is no need to look further
		if(bestVertexPair.second == 0) {
			std::cout << "no weight gain. Get buffed, bro" << std::endl;
			break;
		}
		// add most weight Vertex to Solution
		greedySolution.insert(bestVertexPair.first);
		currentWeight += bestVertexPair.second;

	}
	// Vertices with degree of 1 in the greedy subgraph can be made longer without compromising the other Nodes
	// construct the subgraph with all nodes from/to solution nodes in a Segment
	undirectedAdjGraph greedySubGraph = getSolutionSubgraph(greedySolution);

	std::vector<Vertex> deg1Vertices = getDeg1Vertices(greedySubGraph);
	// try to pin other Vertex in the same segment and augment length in original graph
	for(Vertex v : deg1Vertices) {
		//get other Vertex in the same segmentVertex otherVertex;
		auto segmentVertexPair = getOtherSolVertexInSegment(greedySolution, v);
		if (segmentVertexPair.second != Vertex() - 1) {
			Vertex maxLengthVector = getMaxLengthVertex(segmentVertexPair.first, segmentVertexPair.second);
			//augment segment
			if (greedySolution.find(maxLengthVector) != greedySolution.end()) {
				//std::cout << "something went wrong and it chose another solution vertex" << std::endl;
			} else {
				greedySolution.erase(v);
				greedySolution.insert(maxLengthVector);
			}
		}
	}

	return computeSolution(greedySolution);
}

Vertex Capturing2D::getMaxLengthVertex(ulong segmentIndex, Vertex otherVertex) const {
	Vertex maxLengthVector = Vertex() - 1;
	NT maxLength = 0;
	Vector2D pinnedVector = vectorMapping.at(otherVertex);
	for (Vector2D vec : vectorFamilies[segmentIndex]) {

				NT length = GeoUtil::getLengthVector2(pinnedVector, vec);
				if (length > maxLength) {
					maxLength = length;
					maxLengthVector = vertexMapping.at(vec);
				}
			}
	return maxLengthVector;
}

std::pair<ulong, Vertex> Capturing2D::getOtherSolVertexInSegment(const std::set<Vertex> &greedySolution, Vertex v) const {
	Vertex otherVertex= Vertex() - 1;
	ulong segmentIndex= 0;
	NT lengthToOV = 0;
	for (unsigned long segIndex : vertexInFamiliesMapping.at(v)) {
		for (Vertex oV : vertexFamilies[segIndex]) {
			if (greedySolution.find(oV) != greedySolution.end()) {
				NT innerLength = GeoUtil::getLengthVector2(vectorMapping.at(v), vectorMapping.at(oV));
				if (innerLength > lengthToOV) {
					segmentIndex = segIndex;
					otherVertex = oV;
					lengthToOV = innerLength;
				}
			}
		}
	}
	return std::make_pair(segmentIndex,otherVertex);
}

std::pair<Vertex, NT>
Capturing2D::getBestVertex(std::set<Vertex> &greedySolution, const std::map<Vector2D, Vertex> &vertexMappingSub,
						   const std::map<Vertex, Vector2D> &vectorMappingSub,
						   const std::vector<std::vector<Vector2D>> &vectorFamiliesSub) {
	NT mostWeightGain= 0;
	Vertex mostWeightVertex= Vertex() - 1;
	std::set<Vertex> possibleVertices = getPossibleVertices(greedySolution, vertexMappingSub, vectorMappingSub,
															vectorFamiliesSub);

	// search for most weight for all possible Vertices//std::map<unsigned long,std::pair<Vertex,Vertex>> newSegmentIndexesMap;
	for(Vertex v : possibleVertices) {
		NT weight = 0;
		//get all segment indexes that are adjacent to vertex v
		const std::set<unsigned long> &segmentIndexes = vectorSegmentMapping[vectorMappingSub.at(v)];
		//std::map<unsigned long,std::pair<Vertex,Vertex>> possibleSegmentIndexesMap;
		for (unsigned long segIndex : segmentIndexes) {
			//when adjacent segment is already in solution there will be a weight gain
			NT innerMaxWeight = 0;
			NT existingSegmentWeight = 0;
			auto familyIndexPair = getFamilyIndexPair(segIndex, greedySolution);
			if (familyIndexPair.second > familyIndexPair.first) {
				existingSegmentWeight = getWeightInFamilyWithIndex(segIndex, familyIndexPair.first,
																   familyIndexPair.second);
			}
			for (auto vec : vectorFamiliesSub[segIndex]) {
				if (greedySolution.find(vertexMappingSub.at(vec)) != greedySolution.end()) {
					NT length = GeoUtil::getLengthVector2(vec, vectorMappingSub.at(v));
					if ((length - existingSegmentWeight) > innerMaxWeight) {
						innerMaxWeight = length;
					}
				}
			}
			weight += innerMaxWeight;
		}
		if (weight > mostWeightGain) {
			mostWeightGain = weight;
			mostWeightVertex = v;
		}
	}
	return std::make_pair(mostWeightVertex,mostWeightGain);
}

std::set<Vertex> Capturing2D::getPossibleVertices(const std::set<Vertex> &greedySolution,
												  const std::map<Vector2D, Vertex> &vertexMappingSub,
												  const std::map<Vertex, Vector2D> &vectorMappingSub,
												  const std::vector<std::vector<Vector2D>> &vectorFamiliesSub) const {
	std::set<Vertex> possibleVertices;
	for(Vertex v : greedySolution) {
			for(Vertex segIndex : vectorSegmentMapping.at(vectorMappingSub.at(v))) {
				for (auto vec : vectorFamiliesSub.at(segIndex)) {
					possibleVertices.insert(vertexMappingSub.at(vec));
				}
			}
		}
	//delete already used vertices afterwards
	for(Vertex v : greedySolution) {
			possibleVertices.erase(v);
		}
	return possibleVertices;
}

void Capturing2D::AddToSegmentIndexesToUsedVertices(const std::map<Vector2D, Vertex> &vertexMappingSub,
													const Vector2D &vector,
													std::map<Vertex, std::set<Vertex>> *segmentIndexesToUsedVertices) const {
	for(unsigned long segIndex : vectorSegmentMapping.at(vector)) {
		if(segmentIndexesToUsedVertices->find(segIndex) == segmentIndexesToUsedVertices->end()) {
			segmentIndexesToUsedVertices->insert(
					std::pair<unsigned int, std::set<Vertex>>(segIndex, std::set<Vertex>()));
		}
		(*segmentIndexesToUsedVertices)[segIndex].insert(vertexMappingSub.at(vector));
	}
}

std::vector<Vertex> Capturing2D::getDeg1Vertices(const undirectedAdjGraph &greedySubGraph) const {
	std::vector<Vertex> deg1Vertices;
	//look for Vertices with degree of 1
	for(Vertex v : greedySubGraph.getVertices()) {
		if(greedySubGraph.getEdges(v).size() <= 1) {
			deg1Vertices.push_back(v);
		}
	}
	return deg1Vertices;
}

undirectedAdjGraph Capturing2D::getSolutionSubgraph(const std::set<Vertex> &greedySolution) {
	std::set<Vertex> allVertexSolution;
	for (unsigned long segIndex = 0; segIndex < segments.size(); ++segIndex) {
		auto indexes = getFamilyIndexPair(segIndex, greedySolution);
		unsigned long j = indexes.first;
		unsigned long k = indexes.second;
		//found a segment track
		if(k > j) {
			for (unsigned long l = j; l <= k ; ++l) {
				// add every vertex to the new subgraph
				allVertexSolution.insert(vertexFamilies[segIndex][l]);
			}
		}
	}
	return graph.createSubgraph(allVertexSolution);
}

std::pair<std::pair<Vector2D, Vector2D>, NT>
Capturing2D::getHeaviestFamily(const std::vector<std::vector<Vector2D>> &vectorFamiliesSub) const {
	NT longestEdgeLength= 0;
	Vector2D longestVectorPair[2];
	for (unsigned long j = 0; j < vectorFamiliesSub.size(); ++j) {
		if(vectorFamiliesSub[j].size() > 1) {
			NT edgeLength = GeoUtil::getLengthVector2(vectorFamiliesSub[j][0],
														   vectorFamiliesSub[j][vectorFamiliesSub[j].size() - 1]);
			if (edgeLength > longestEdgeLength) {
				longestEdgeLength = edgeLength;

				longestVectorPair[0] = vectorFamiliesSub[j][0];
				longestVectorPair[1] = vectorFamiliesSub[j][vectorFamiliesSub[j].size() - 1];
			}
		}
	}
	return std::make_pair(std::make_pair(longestVectorPair[0],longestVectorPair[1]),longestEdgeLength);
}

undirectedAdjGraph Capturing2D::computeConnectedSubgraph(std::map<Vector2D, Vertex> *vertexMappingSub,
														 std::map<Vertex, Vector2D> *vectorMappingSub,
														 std::vector<std::vector<Vector2D>> *vectorFamiliesSub) const {
	undirectedAdjGraph subGraph = graph;
	std::set<Vertex> connectedVertices;
	while(connectedVertices.size() != subGraph.getSize()) {
		connectedVertices.clear();
		vectorMappingSub->clear();
		vertexMappingSub->clear();
		for (Vertex v : subGraph.getVertices()) {
			if (subGraph.getEdges(v).size() > 1) {
				connectedVertices.insert(v);
				vertexMappingSub->insert(std::pair<Vector2D, Vertex>(vectorMapping.at(v), v));
				vectorMappingSub->insert(std::pair<Vertex, Vector2D>(v, vectorMapping.at(v)));
			}
		}
		subGraph = subGraph.createSubgraph(connectedVertices);
	}
	//translate the points from original graph to subgraph
	for (unsigned int i = 0; i < vectorFamilies.size(); ++i) {
		vectorFamiliesSub->push_back(std::vector<Vector2D>());
		for (Vector2D vec : vectorFamilies[i]) {
			if(connectedVertices.find(vertexMapping.at(vec)) != connectedVertices.end()) {
				(*vectorFamiliesSub)[i].push_back(vec);
			}
		}
	}
	return subGraph;
}

void Capturing2D::setBounds(uint xBound, uint yBound) {
	this->xBound = xBound;
	this->yBound = yBound;
}

void Capturing2D::computeIntersectionPoints(Segments *segments, std::set<Point_2> *points,
											std::vector<std::set<Point_2>> *familyPoints,
											std::map<Point_2, std::set<Vertex>> *vectorVertexMap) {
	// Create the corresponding vector of bounding boxes
	std::vector<Box> boxes;
	for ( Iterator i = segments->begin(); i != segments->end(); ++i) {
		Box b = Box(i->bbox(), i);
		boxes.push_back(b);
	}


	// Create the corresponding vector of pointers to bounding boxes
	std::vector<Box *> ptr;
	for ( std::vector<Box>::iterator i = boxes.begin(); i != boxes.end(); ++i)
		ptr.push_back(&*i);


	// Run the self intersection algorithm with all defaults on the
	// indirect pointers to bounding boxes. Avoids copying the boxes.

	CGAL::box_self_intersection_d( ptr.begin(), ptr.end(), CGAL_CapturingReport(segments,familyPoints,vectorVertexMap,points));

}

void Capturing2D::convertToCGALSegments(const std::vector<Segment<NT>> &segs, Segments *cgal_segments) {
	cgal_segments->clear();
	for (auto s : segs) {
		Point_2 cgal_start(s.GetStartVector().x,s.GetStartVector().y);
		Point_2 cgal_end(s.GetEndVector().x,s.GetEndVector().y);


		cgal_segments->push_back(Segment_2(cgal_start,cgal_end));
	}
}

void Capturing2D::convertCGALToGraph(const Segments &segments, std::set<Point_2> *points,
									 std::vector<std::set<Point_2>> *familyPoints,
									 std::map<Point_2, std::set<Vertex>> *vectorVertexMap, std::set<Vector2D> *vecSet) {

	for (uint i = 0; i < segments.size() ; ++i) {
		if(points->find(segments[i].source()) == points->end()) {
			points->insert(segments[i].source());
		}
		if(points->find(segments[i].target()) == points->end()) {
			points->insert(segments[i].target());
		}
	}

	//translate points to vectorFamilies
	if(vectorFamilies.size() != familyPoints->size()) {
		std::cout << "vectorFamilies and reportFamilies does not have the same size: "
				  << vectorFamilies.size() << " to " << familyPoints->size() << std::endl;
	}
	//vectorFamilies.clear();
	for (uint j = 0; j < familyPoints->size(); ++j) {
		//vectorFamilies.push_back(std::vector<Vector2D>());
		//auto start = segments[j].GetStartVector();
		//auto end = segments[j].GetEndVector();
		//vectorFamilies[j].push_back(start);
		auto segStart = segments[j].source();
		auto segEnd = segments[j].target();
		if((*familyPoints)[j].find(segStart) == (*familyPoints)[j].end()) {
			(*familyPoints)[j].insert(segStart);
		}
		if((*familyPoints)[j].find(segEnd) == (*familyPoints)[j].end()) {
			(*familyPoints)[j].insert(segEnd);
		}
		for (Point_2 point_2 : (*familyPoints)[j]) {
			Vector2D vec(CGAL::to_double(point_2.x()), CGAL::to_double(point_2.y()));
			vectorFamilies[j].push_back(vec);
		}
		//vectorFamilies[j].push_back(end);
	}
	for (auto point_pair : *vectorVertexMap) {
		std::pair<Vector2D , std::set<Vertex>> vec_pair(Vector2D(CGAL::to_double(point_pair.first.x()),CGAL::to_double(point_pair.first.y()))
														,std::set<Vertex>());
		vec_pair.second.insert(point_pair.second.begin(),point_pair.second.end());
		vectorSegmentMapping.insert(vec_pair);
	}
	vecSet->clear();
	for(auto point : *points) {
		vecSet->insert(Vector2D(CGAL::to_double(point.x()),CGAL::to_double(point.y())));
	}
}

template<class T>
inline bool vectorContains(std::vector<T> *vector, const T element) {
	return std::find(vector->begin(), vector->end(), element) != vector->end();
}

std::set<Vector2D> * Capturing2D::calculateSegmentIntersection() {
	std::set<Vector2D>* vecSet = new std::set<Vector2D>();
    for (unsigned int i = 0; i < segments.size(); ++i) {
        Segment<NT> seg1 = segments[i];

        //std::vector<Vector2D> *seg1Points = &vectorFamilies[i];
		vecSet->insert((Vector2D &&) seg1.GetStartVector());
		vecSet->insert((Vector2D &&) seg1.GetEndVector());
    	//auto k = seg1Points->begin();
		if(!vectorContains(&vectorFamilies[i],seg1.GetStartVector())) {
			vectorFamilies[i].push_back((Vector2D &&)seg1.GetStartVector());

		} else {
			std::cout << "found already inserted Vector " << seg1.GetStartVector().x << "," << seg1.GetStartVector().y << std::endl;
		}
		if(!vectorContains(&vectorFamilies[i],seg1.GetEndVector())) {
			vectorFamilies[i].push_back((Vector2D &&)seg1.GetEndVector());
		}else {
			std::cout << "found already inserted Vector " << seg1.GetEndVector().x << "," << seg1.GetEndVector().y << std::endl;
		}

		if(vectorSegmentMapping.find(seg1.GetStartVector()) == vectorSegmentMapping.end()) {
			vectorSegmentMapping.insert(
					std::pair<Vector2D, std::set<Vertex>>(seg1.GetStartVector(), std::set<Vertex>()));
		}
		if(vectorSegmentMapping.find(seg1.GetEndVector()) == vectorSegmentMapping.end()) {
			vectorSegmentMapping.insert(
					std::pair<Vector2D, std::set<Vertex>>(seg1.GetEndVector(), std::set<Vertex>()));
		}
		vectorSegmentMapping[seg1.GetStartVector()].insert(i);
		vectorSegmentMapping[seg1.GetEndVector()].insert(i);

		//check all other segments for intersection points
        for (unsigned int j = i+1; j < segments.size(); ++j) {
            Segment<NT> seg2 = segments[j];
            std::pair<Vector2D, bool> intersectPair = GeoUtil::getIntersectionPoint(seg1,seg2);
            if(intersectPair.second) {
                //std::vector<Vector2D> *seg2Points = &vectorFamilies[j];
				if(!vectorContains(&vectorFamilies[i],intersectPair.first)) {
					vectorFamilies[i].push_back(intersectPair.first);
				}else {
					std::cout << "found already inserted Vector " << intersectPair.first.x << "," << intersectPair.first.y << std::endl;
				}
				if(!vectorContains(&vectorFamilies[j],intersectPair.first)) {
					vectorFamilies[j].push_back(intersectPair.first);
				} else {
					std::cout << "found already inserted Vector " << intersectPair.first.x << "," << intersectPair.first.y << std::endl;
				}
				vecSet->insert(intersectPair.first);
				segmentGraph.setEdge(i,j, 1);
				if(vectorSegmentMapping.find(intersectPair.first) == vectorSegmentMapping.end()) {
					vectorSegmentMapping.insert(
							std::pair<Vector2D, std::set<Vertex>>(intersectPair.first, std::set<Vertex>()));
				}
				vectorSegmentMapping[intersectPair.first].insert(i);
				vectorSegmentMapping[intersectPair.first].insert(j);
            }
        }
    }
	return vecSet;
}

void Capturing2D::translateGraph(const std::set<Vector2D> &pSet) {
	vertexMapping.clear();
	vectorMapping.clear();
	graph.clear();

	for(auto vec : pSet) {
		Vertex v = graph.addVertex();
		vertexMapping.insert(std::pair<Vector2D,Vertex>(vec,v));
		vectorMapping.insert(std::pair<Vertex,Vector2D>(v,vec));
	}
	// translate vector families
	vertexFamilies = std::vector<std::vector<Vertex>>(vectorFamilies.size());
	for (uint j = 0; j < vectorFamilies.size(); ++j) {
		std::sort(vectorFamilies[j].begin(),vectorFamilies[j].end());
		Vertex v = vertexMapping[vectorFamilies[j][0]];
		vertexFamilies[j].push_back(v);
		for (uint k = 1; k < vectorFamilies[j].size(); ++k) {
			Vertex v2 = vertexMapping[vectorFamilies[j][k]];
			vertexFamilies[j].push_back(v2);
			NT weight = GeoUtil::getLengthVector2(vectorFamilies[j][k-1],vectorFamilies[j][k]);
			auto edgeResult = graph.setEdge(v, v2, CGAL::to_double(weight));
			//std::cout << "Vector familie: " <<  vectorFamilies[j][k].x << "," << vectorFamilies[j][k].y << std::endl;
			//std::cout << "Family: "<< j << " Edge: " << v << "," << v2 << ": " << weight << std::endl;
			//std::cout << std::fixed << std::setprecision(20) << "Family: "<< j << " Edge: " << vectorFamilies[j][k-1] << ";" << vectorFamilies[j][k-1] << std::endl;
			if(!edgeResult) {
				std::cout << "Fehler beim Eintragen in die Adjazenzliste bei " << v << "," << v2 << std::endl;
			}
			v = v2;
		}
	}
	/*
	std::map<Segment, std::vector<Vertex>> segmentToVertexMap;
	for (unsigned int k = 0; k < segments.size(); ++k) {
		// first sort the small lists
		Segment seg = segments[k];
		std::vector<Vector2D> *segVecs = &vectorFamilies[k];
		std::sort(segVecs->begin(),segVecs->end());

		Vertex vertice1;
		auto iterator = vertexMapping.find(segVecs->at(0));
		if(iterator == vertexMapping.end()) {
			vertice1 = graph.addVertex();
			vertexMapping.insert(std::pair<Vector2D, Vertex>(segVecs->at(0), vertice1));
			vectorMapping.insert(std::pair<Vertex, Vector2D>(vertice1, segVecs->at(0)));
		} else {
			vertice1 = iterator->second;
		}
		segmentToVertexMap.insert(std::pair<Segment, std::vector<Vertex>>(seg,std::vector<Vertex>()));
		segmentToVertexMap[seg].push_back(vertice1);
		vertexFamilies[k].push_back(vertice1);

		for (unsigned int i = 1; i < segVecs->size(); ++i) {
			Vertex vertice2;
			auto mappingIterator = vertexMapping.find(segVecs->at(i));
			if(mappingIterator == vertexMapping.end()) {
				vertice2 = graph.addVertex();
				vertexMapping.insert(std::pair<Vector2D, Vertex>(segVecs->at(i), vertice2));
				vectorMapping.insert(std::pair<Vertex, Vector2D>(vertice2, segVecs->at(i)));
			} else {
				vertice2 = mappingIterator->second;
			}
			segmentToVertexMap[seg].push_back(vertice2);
			vertexFamilies[k].push_back(vertice2);
			auto edgeResult = graph.setEdge(vertice1, vertice2, GeoUtil::getLengthVector2D(segVecs->at(i-1),segVecs->at(i)));
			if(!edgeResult) {
				std::cout << "Fehler beim Eintragen in die Adjazenzliste bei " << vertice1 << "," << vertice2 << std::endl;
			}
			vertice1 = vertice2;
		}
	}
	*/


	// map vector segments to vertex family indexes
	for (auto l = vectorSegmentMapping.begin(); l != vectorSegmentMapping.end(); ++l) {
		vertexInFamiliesMapping.insert(std::pair<Vertex,std::set<Vertex>>(vertexMapping[l->first],std::set<Vertex>()));
		vertexInFamiliesMapping[vertexMapping[l->first]].insert(l->second.begin(),l->second.end());
	}
}

std::pair<std::pair<Vector2D, Vector2D>, NT>
Capturing2D::getFamilyMaxLength(ulong segment, std::vector<std::vector<Vector2D>> vectorFamiliesSub) {
	NT longestEdgeLength = 0;
	Vector2D longestVectorPair[2];
	if (vectorFamiliesSub[segment].size() > 1) {
		NT edgeLength = GeoUtil::getLengthVector2(vectorFamiliesSub[segment][0],
													   vectorFamiliesSub[segment][vectorFamiliesSub[segment].size() -
																				  1]);
		longestEdgeLength = edgeLength;

		longestVectorPair[0] = vectorFamiliesSub[segment][0];
		longestVectorPair[1] = vectorFamiliesSub[segment][vectorFamiliesSub[segment].size() - 1];

		return std::make_pair(std::make_pair(longestVectorPair[0], longestVectorPair[1]), longestEdgeLength);
	} else {
		return std::make_pair(std::make_pair(Vector2D(0,0), Vector2D(0,0)), longestEdgeLength);
	}

}





