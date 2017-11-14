//
// Created by gharg on 30.08.17.
//

#ifndef PROJECT_SVGIMPORT_H
#define PROJECT_SVGIMPORT_H

#include <vector>
#include "../Utilities/GeoUtil.h"
#include "../Capturing/Header/Capturing2D.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>

struct Polyline {
	std::vector<Vector2D> points;
};

struct Line {
	Vector2D vec1;
	Vector2D vec2;
};

class Imports {

//	read_xml(is, pt);
public:
	static Capturing2D importSVG(std::istream & istream);

	static std::vector<Vector2D> importCSV(std::ifstream &ifstream, const char *lineseparator,
										   const char *columnseparator);
};


#endif //PROJECT_SVGIMPORT_H
