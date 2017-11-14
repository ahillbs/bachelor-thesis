//
// Created by gharg on 19.09.17.
//

#ifndef PROJECT_GENERATORWRAPPER_H
#define PROJECT_GENERATORWRAPPER_H

#include "../Header/Capturing2D.h"
#include "../../MacroTestingFramework/TestInstance.hpp"
#include "../Generator/Header/Cap2DGenerator.h"

static std::vector<Vector2D> transformTestInstanceToVector2D(const std::vector<std::vector<NT>> &points,ulong amount) {
	std::vector<Vector2D> vecs;
	for (uint i = 0; i < amount; ++i) {
		//std::cout << "points: " << points[0][i] << " " << points[1][i] << std::endl;
		Vector2D v(points[0][i],points[1][i]);
		vecs.push_back(v);
	}
	return vecs;
}

static TestInstance<NT> transformVector2DToTestInstance(const std::vector<Vector2D> &points) {
	std::vector<NT> vecs;
	for (uint i = 0; i < points.size(); ++i) {
		//std::cout << "points: " << points[0][i] << " " << points[1][i] << std::endl;
		vecs.push_back(points[i].x);
		vecs.push_back(points[i].y);
	}
	return TestInstance<NT>("randomPoints", &vecs[0], &vecs[vecs.size()], 2);
}

class RandomSharedEndPointsInstanceGenerator {
	ulong segmentAmount = 100;
	uint maxK = 5;
	bool isRelativeSegmentAmount = false;

public:

	ulong getSegmentAmount() const {
		return segmentAmount;
	}

	void setSegmentAmount(ulong segmentAmount) {
		this->segmentAmount = segmentAmount;
	}

	uint getMaxK() const {
		return maxK;
	}

	void setMaxK(uint maxK) {
		this->maxK = maxK;
	}

	bool isIsRelativeSegmentAmount() const {
		return isRelativeSegmentAmount;
	}

	void setIsRelativeSegmentAmount(bool isRelativeSegmentAmount) {
		this->isRelativeSegmentAmount = isRelativeSegmentAmount;
	}

	Capturing2D generateInstance(const std::vector<std::vector<NT>> &points,ulong amount,ulong segments ,ulong seedStart) {
		auto vecs = transformTestInstanceToVector2D(points,amount);
		ulong totalSegmentAmount = calcSegmentAmount(vecs.size(), segments);
		Capturing2D cap2D = Cap2DGenerator::generateRandomSharedEndPointsInstance(vecs,totalSegmentAmount,seedStart,maxK);
		cap2D.convertSegments2Graph();
		return cap2D;
	}

	ulong calcSegmentAmount(ulong pointAmount, ulong segments) const {
		ulong totalSegmentAmount = segments;
		ulong maxSegments = (pointAmount * (pointAmount-1)) / 2;
		if(isRelativeSegmentAmount) {
			if(segments >= 100) {
				totalSegmentAmount = maxSegments;
			} else {
				double percent = maxSegments * ((double) segments / (double)100);
				totalSegmentAmount = (ulong)(percent);
			}
		} else {
			if(segments > maxSegments) {
				totalSegmentAmount = maxSegments;
			}
		}
		return totalSegmentAmount;
	}


	std::string printSettings() {
		std::string result = "maxK = ";
		result += std::to_string(maxK);
		result += ", ";
		if(isRelativeSegmentAmount) {
			result += "segment percentage ";
			result += std::to_string(segmentAmount);
			result += "%";
		} else {
			result += "segments: ";
			result += std::to_string(segmentAmount);
		}
		return result;
	}

};

class ProbabilisticInstanceGenerator {
	double probability = 20;
	ulong maxK = 5;

public:


	ulong getMaxK() const {
		return maxK;
	}

	void setMaxK(ulong maxK) {
		this->maxK = maxK;
	}

	double getSegmentAmount() const {
		return probability;
	}
	void setSegmentAmount(double prob) {
		this->probability = prob;
	}
	double getProbability() const {
		return probability;
	}

	void setProbability(double probability) {
		this->probability = probability;
	}

	Capturing2D generateInstance(const std::vector<std::vector<NT>> &points,ulong amount,ulong segments ,ulong seedStart) {
		auto vecs = transformTestInstanceToVector2D(points,amount);

		Capturing2D cap2D = Cap2DGenerator::generateProbabilisticInstance(vecs,segments,seedStart,maxK);
		cap2D.convertSegments2Graph();
		return cap2D;
	}

	std::string printSettings() {
		std::string result = "maxK = ";
		result += std::to_string(maxK);
		result += ", ";
		result += "probability ";
		result += std::to_string(probability);
		result += "%";

		return result;
	}
};

class RandomAllNodesInstanceGenerator {
	ulong segmentAmount = 100;
	uint maxK = 5;
	bool isRelativeSegmentAmount = false;

public:

	ulong getSegmentAmount() const {
		return segmentAmount;
	}

	void setSegmentAmount(ulong segmentAmount) {
		this->segmentAmount = segmentAmount;
	}

	uint getMaxK() const {
		return maxK;
	}

	void setMaxK(uint maxK) {
		this->maxK = maxK;
	}

	bool isIsRelativeSegmentAmount() const {
		return isRelativeSegmentAmount;
	}

	void setIsRelativeSegmentAmount(bool isRelativeSegmentAmount) {
		this->isRelativeSegmentAmount = isRelativeSegmentAmount;
	}

	Capturing2D generateInstance(const std::vector<std::vector<NT>> &points,ulong amount,ulong segments ,ulong seedStart) {
		auto vecs = transformTestInstanceToVector2D(points,amount);
		ulong totalSegmentAmount = calcSegmentAmount(vecs.size(), segments);
		Capturing2D cap2D = Cap2DGenerator::generateRandomAllNodesInstance(vecs,totalSegmentAmount,seedStart,maxK);
		cap2D.convertSegments2Graph();
		return cap2D;
	}

	ulong calcSegmentAmount(ulong pointAmount, ulong segments) const {
		ulong totalSegmentAmount = segments;
		ulong maxSegments = (pointAmount * (pointAmount-1)) / 2;
		if(isRelativeSegmentAmount) {
			if(segments >= 100) {
				totalSegmentAmount = maxSegments;
			} else {
				double percent = maxSegments * ((double) segments / (double)100);
				totalSegmentAmount = (ulong)(percent);
			}
		} else {
			if(segments > maxSegments) {
				totalSegmentAmount = maxSegments;
			}
		}
		return totalSegmentAmount;
	}
	std::string printSettings() {
		std::string result = "maxK = ";
		result += std::to_string(maxK);
		result += ", ";
		if(isRelativeSegmentAmount) {
			result += "segment percentage ";
			result += std::to_string(segmentAmount);
			result += "%";
		} else {
			result += "segments: ";
			result += std::to_string(segmentAmount);
		}
		return result;
	}
};

#endif //PROJECT_GENERATORWRAPPER_H
