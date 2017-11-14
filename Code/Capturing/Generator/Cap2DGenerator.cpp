//
// Created by gharg on 09.09.17.
//

#include <random>
#include <functional>
#include "Header/Cap2DGenerator.h"
#include "../../Visualization/Imports.h"


std::vector<Vector2D> Cap2DGenerator::getRandomPoints(uint amount, uint xBound, uint yBound, ulong seed) {

	std::mt19937_64 random(seed);
	std::mt19937_64 randomY(seed+1);

	auto dice_randX = std::bind(std::uniform_int_distribution<int>(0,xBound), random);
	auto dice_randY = std::bind(std::uniform_int_distribution<int>(0,yBound), randomY);
	std::vector<Vector2D> vectors;
	for (unsigned long i = 0; i < amount; ++i) {
		int coordinates[2];
		coordinates[0] = dice_randX();
		coordinates[1] = dice_randY();
		vectors.push_back(Vector2D(coordinates[0],coordinates[1]));
	}
	return vectors;
}

Capturing2D
Cap2DGenerator::generateRandomSharedEndPointsInstance(const std::vector<Vector2D> &points, ulong segmentAmount,
													  ulong seed, uint k) {
	std::vector<Segment<NT>> segments;
	std::mt19937_64 randomS(seed+2);
	NT maxX = 0,maxY = 0;
	unsigned long maxEdges = points.size()*(points.size()-1) / 2;
	/*std::vector<Segment> allSegments;
	for (unsigned long j = 0; j < maxEdges; ++j) {
		for (unsigned long i = j+1; i < maxEdges; ++i) {
			allSegments.push_back(Segment(vectors[i], vectors[j]));
		}
	}*/

	if(segmentAmount < maxEdges) {
		auto dice_randS = std::bind(std::uniform_int_distribution<unsigned long>(0, points.size() - 1), randomS);

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
				if(maxX < points[vecIndex1].x || maxX < points[vecIndex2].x) {
					maxX = points[vecIndex1].x > points[vecIndex2].x ?  points[vecIndex1].x : points[vecIndex2].x;
				}
				if(maxY < points[vecIndex1].y || maxX < points[vecIndex2].y) {
					maxY = points[vecIndex1].y > points[vecIndex2].y ?  points[vecIndex1].y : points[vecIndex2].y;
				}
				longPair pair2 = longPair(vecIndex2,vecIndex1);
				indexPairs.insert(pair);
				indexPairs.insert(pair2);

				segments.push_back(Segment<NT>(points[vecIndex1], points[vecIndex2]));
			} else {

			}
		}
	} else {
		std::cout << "Amount of Segments exceed maximum Edges for " << points.size() << " points! " << std::endl;
		std::cout << "Complete graph will be build instead" << std::endl;
		for (unsigned long j = 0; j < maxEdges; ++j) {
			if(maxX < points[j].x ) {
				maxX = points[j].x;
			}
			if(maxY < points[j].y ) {
				maxY = points[j].y;
			}
			for (unsigned long i = j+1; i < maxEdges; ++i) {
				segments.push_back(Segment<NT>(points[i], points[j]));
			}
		}
	}
	return Capturing2D(k, segments, (uint) (CGAL::to_double(maxX) * 1.05), (uint) (CGAL::to_double(maxY) * 1.05));
}

Capturing2D Cap2DGenerator::generateRandomSegmentInstance(uint amount, uint xBound, uint yBound, ulong seed, uint k = 5) {

	std::mt19937_64 random(seed);
	std::mt19937_64 randomY(seed+1);

	std::vector<Segment<NT>> segments;


	auto dice_randX = std::bind(std::uniform_int_distribution<int>(0,xBound), random);
	auto dice_randY = std::bind(std::uniform_int_distribution<int>(0,yBound), randomY);
	for (uint i = 0; i < amount; ++i) {
		int coordinates[4];
		for (int j = 0; j < 4; ++j) {

			if(j % 2 == 0) {
				coordinates[j] = dice_randX();
			} else {
				coordinates[j] = dice_randY();
			}
		}
		segments.push_back(Segment<NT>(coordinates[0],coordinates[1],coordinates[2],coordinates[3]));
	}
	return Capturing2D(k, segments, xBound, yBound);
}

Capturing2D
Cap2DGenerator::generateProbabilisticInstance(const std::vector<Vector2D> &points, double probability, ulong seed,
											  uint k) {


	std::mt19937_64 random(seed);

	std::vector<Segment<NT>> segments;
	auto dice_rand = std::bind(std::uniform_int_distribution<int>(0,100), random);

	NT maxX = 0;
	NT maxY = 0;
	std::set<std::pair<Vector2D,Vector2D>> pointPair;
	for (uint i = 0; i < points.size(); ++i) {
		if(maxX < points[i].x) {
			maxX = points[i].x;
		}
		if(maxY < points[i].y) {
			maxY = points[i].y;
		}
		for (uint j = i+1; j < points.size(); ++j) {
			if(dice_rand() <= probability) {
				segments.push_back(Segment<NT>(points[i],points[j]));
				auto inserted = pointPair.insert(std::pair<Vector2D,Vector2D>(points[i],points[j]));
				if(!inserted.second) {
					std::cout << points[i] << " and " << points[j] << "was already inserted" << std::endl;
				}
				auto insertedR = pointPair.insert(std::pair<Vector2D,Vector2D>(points[j],points[i]));
				if(!insertedR.second) {
					std::cout << points[i] << " and " << points[j] << "was already inserted" << std::endl;
				}
			}
		}
	}

	return Capturing2D(k, segments, (uint) (CGAL::to_double(maxX) * 1.05), (uint) (CGAL::to_double(maxY) * 1.05));
}


std::vector<Point_22> Cap2DGenerator::calcIntersectionPoints(Segments2& firstSegments,Segments2& intersectionSegments) {
	std::vector<Box2> boxes;

	for ( Iterator2 i = firstSegments.begin(); i != firstSegments.end(); ++i) {
		Box2 b = Box2(i->bbox(), i);
		boxes.push_back(b);
	}
	// Create the corresponding vector of pointers to bounding boxes
	std::vector<Box2 *> ptr;
	for(auto& i: boxes){ptr.push_back(&i);}
	//for ( std::vector<Box>::iterator i = boxes.begin(); i != boxes.end(); ++i)
	//	ptr.push_back(&*i);
	std::vector<Box2> boxes2;
	for ( Iterator2 i = intersectionSegments.begin(); i != intersectionSegments.end(); ++i) {
		Box2 b = Box2(i->bbox(), i);
		boxes2.push_back(b);
	}
	std::vector<Box2 *> ptr2;
	for (std::vector<Box2>::iterator i = boxes2.begin(); i != boxes2.end(); ++i) {
		ptr2.push_back(&(*i));
	}
	// Run the self intersection algorithm with all defaults on the
	// indirect pointers to bounding boxes. Avoids copying the boxes.
	std::vector<Point_22> points;
	std::ptrdiff_t cutoff = 50;

	CGAL::box_intersection_d( ptr.begin(), ptr.end(),ptr2.begin(),ptr2.end(), CGAL_Report(points),cutoff);
	return points;
}

template<class T>
inline bool vectorContains(std::vector<T> *vector, const T element) {
	return std::find(vector->begin(), vector->end(), element) != vector->end();
}

Capturing2D
Cap2DGenerator::generateRandomAllNodesInstance(std::vector<Vector2D> points, ulong segmentAmount, ulong seed,
											   uint k) {
	Segments2 segments;
	//translate points to CGAL points
	std::set<Point_22> points_2;
	Point_22 bBoxPointLower(std::numeric_limits<double>::max(),std::numeric_limits<double>::max());
	Point_22 bBoxPointUpper(std::numeric_limits<double>::min(),std::numeric_limits<double>::min());
	//Arrangement_2 arr;
	//Walk_pl pl(arr);
	for (uint l = 0; l < points.size(); ++l) {
		Point_22 p = Point_22(CGAL::to_double(points[l].x),CGAL::to_double(points[l].y));
		points_2.insert(p);
		if(bBoxPointLower.x() > p.x()) {
			bBoxPointLower = Point_22(p.x(),bBoxPointLower.y());
		}
		if(bBoxPointLower.y() > p.y()) {
			bBoxPointLower = Point_22(bBoxPointLower.x(),p.y());
		}
		if(bBoxPointUpper.x() < p.x()) {
			bBoxPointUpper = Point_22(p.x(),bBoxPointUpper.y());
		}
		if(bBoxPointUpper.y() < p.y()) {
			bBoxPointUpper = Point_22(bBoxPointUpper.x(),p.y());
		}

		//CGAL::insert_point(arr,p);
	}
	//std::cout << bBoxPointLower << " ... " << bBoxPointUpper << std::endl;
	//std::cout << CGAL::to_double(bBoxPointLower.x()) << " :.: " << CGAL::to_double(bBoxPointLower.y()) << std::endl;
	//std::cout << CGAL::to_double(bBoxPointUpper.x()) << " .:. " << CGAL::to_double(bBoxPointUpper.y()) << std::endl;


	NT2 maxX = 10, maxY = 10;
	//unsigned long maxEdges = points.size() * (points.size() - 1) / 2;
	//std::cout << segmentAmount << " segments will be generated" << std::endl;
	using longPair = std::pair<unsigned long, unsigned long>;
	std::set<longPair> indexPairs;
	for (unsigned long j = 0; j < segmentAmount; ++j) {
		//std::cout << "Segment nr.: " << j << std::endl;
		std::mt19937_64 randomS(seed + 2 + points_2.size() + j);
		auto dice_randS = std::bind(std::uniform_int_distribution<unsigned long>(0, points_2.size() - 1), randomS);
		auto vecIndex1 = dice_randS();
		auto vecIndex2 = dice_randS();
		while (vecIndex1 == vecIndex2) {
			vecIndex2 = dice_randS();
		}
		longPair pair = longPair(vecIndex1, vecIndex2);
		auto it1 = points_2.begin();//arr.vertices_begin();
		auto it2 = points_2.begin();//arr.vertices_begin();
		if (indexPairs.find(pair) == indexPairs.end()) {
			std::advance(it1 , vecIndex1);
			std::advance(it2 , vecIndex2);
			if (maxX < it1->x() || maxX < it2->x()) {
				maxX = it1->x() > it2->x() ? it1->x() : it2->x();
			}
			if (maxY < it1->y() || maxX < it2->y()) {
				maxY = it1->y() > it2->y() ? it1->y() : it2->y();
			}
			longPair pair2 = longPair(vecIndex2, vecIndex1);
			indexPairs.insert(pair);
			indexPairs.insert(pair2);
			Segment_22 loneSegment(*it1, *it2);
			//std::cout << "segment " << loneSegment << std::endl;
			//std::cout << "point 1: " << *it1 << " point 2: " << *it2 << std::endl;
			//Segments2 loneSeg = {Segment_22(it1->point(), it2->point())};
			//std::cout << "calc intersections" << std::endl;
			//CGAL::insert(arr,loneSegment,pl,)
			//std::vector<Point_22> interPoints =
			//		calcIntersectionPoints(segments,loneSeg);
			//CGAL::insert(arr,loneSegment,pl);
			//std::cout << "new segment " << loneSegment << std::endl;
			for (auto l = segments.begin(); l != segments.end(); ++l) {
				//std::cout << "get intersection " << *l << std::endl;
				auto result = CGAL::intersection(*l,loneSegment);
				if(result.is_initialized()) {
					if (const Point_22 *p = boost::get<Point_22>(&*result)) {
						if(p->x() >= bBoxPointLower.x() && p->y() >= bBoxPointLower.y() &&
						   p->x() <= bBoxPointUpper.x() && p->y() <= bBoxPointUpper.y()) {
							points_2.insert(*p);
						}
						//else
							//std::cout << "found vertex out of bounds" << *p << ":" << dx << " " << dy << std::endl;
					} else {
						//const Segment_2 *s = boost::get<Segment_2>(&*result);
						//std::cout << "found Segment: " << *s << std::endl;
					}
				}
			}
			//std::cout << "finished" << std::endl;
			/*for (auto it = interPoints.begin(); it != interPoints.end(); ++it) {
				if(points_2.find(*it) == points_2.end()) {
					std::cout << "added intersection point " << CGAL::to_double(it->x()) << ":" << CGAL::to_double(it->y()) << std::endl;
					points_2.insert(*it);
				}

			}*/
			segments.push_back(loneSegment);
		} else {

		}
	}
	double xBound = CGAL::to_double(maxX);
	double yBound = CGAL::to_double(maxY);
	xBound *= 1.05;
	yBound *= 1.05;
	//std::cout << "finished building segment graph " << segments.size() << std::endl;
	std::vector<Segment<NT>> realSegments;
	for (auto it = segments.begin(); it != segments.end(); ++it) {

		NT2 nx1 = it->source().x();
		NT2 nx2 = it->target().x();
		NT2 ny1 = it->source().y();
		NT2 ny2 = it->target().y();
		/*double x1 = CGAL::to_double(nx1);
		double x2 = CGAL::to_double(nx2);
		double y1 = CGAL::to_double(ny1);
		double y2 = CGAL::to_double(ny2);*/
		//std::cout << nx1 << " " << nx2 << " " << ny1 << " " << ny2 << std::endl;
		//std::cout << x1 << " " << x2 << " " << y1 << " " << y2 << std::endl;
		realSegments.push_back(Segment<NT>(nx1, ny1, nx2, ny2));
	}
	//std::cout << "finished translate segment graph " << realSegments.size() << std::endl;
	return Capturing2D(k, realSegments, (uint) (xBound), (uint) (yBound));
}



