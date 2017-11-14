//
// Created by gharg on 30.08.17.
//

#include "Imports.h"
#include <boost/foreach.hpp>

Capturing2D Imports::importSVG(std::istream &istream) {
	using boost::property_tree::ptree;
	ptree tree;
	try {
		read_xml(istream, tree);
	}catch(const boost::property_tree::ptree_error &e)
	{
		std::cout << e.what() << std::endl;
	}

	std::vector<Segment<NT>> segments;
	double maxX =  std::numeric_limits<double>::min();
	double maxY = std::numeric_limits<double>::min();
	BOOST_FOREACH(auto const& v, tree.get_child("svg") ) {
					ptree::value_type val = ((ptree::value_type)v);
					if(val.first == "line") {
						double x1,x2,y1,y2;
						x1 = val.second.get("<xmlattr>.x1",0);
						x2 = val.second.get("<xmlattr>.x2",0);
						y1 = val.second.get("<xmlattr>.y1",0);
						y2 = val.second.get("<xmlattr>.y2",0);
						Segment<NT> seg(x1,y1,x2,y2);
						segments.push_back(seg);
						if(std::max(x1,x2) > maxX)
							maxX = std::max(x1,x2);

						if(std::max(y1,y2) > maxY)
							maxY = std::max(y1,y2);
					}
					if(val.first == "polyline") {

					}
				}

	return Capturing2D (0, segments, (uint)(maxX *1.05), (uint)(maxY * 1.05));
}
std::vector<Vector2D> Imports::importCSV(std::ifstream &ifstream, const char *lineseparator,
										 const char *columnseparator) {

	std::vector<Vector2D> points;
	if (ifstream) {
		// get length of file:
		ifstream.seekg (0, ifstream.end);
		unsigned long length;
		length = (unsigned long) ifstream.tellg();
		ifstream.seekg (0, ifstream.beg);

		char * buffer = new char [length];

		std::cout << "Reading " << length << " characters... " << std::endl;
		// read data as a block:
		ifstream.read (buffer,length);

		if (ifstream)
			std::cout << "all characters read successfully." << std::endl;
		else
			std::cout << "error: only " << ifstream.gcount() << " could be read" << std::endl;
		ifstream.close();

		// ...buffer contains the entire file...

		std::string csvString(buffer,length);
		// split string into substrings of lines
		unsigned long pos = 0;
		unsigned long endPos = csvString.find(lineseparator,pos);
		if (endPos == std::string::npos) {
			endPos = csvString.size()-1;
		}
		while(pos < csvString.size()) {
			// ignore comment lines
			if(csvString[pos] != '#') {
				unsigned long innerPosMid = csvString.find(columnseparator,pos);
				std::string subStrX = csvString.substr(pos,innerPosMid);
				std::string subStrY = csvString.substr(innerPosMid+1,endPos);
				double x = std::stod(subStrX);
				double y = std::stod(subStrY);
				points.push_back(Vector2D(x,y));
			}


			pos = endPos + 1;
			endPos = csvString.find(lineseparator,pos);
			if (endPos == std::string::npos) {
				endPos = csvString.size()-1;
			}
		}
		delete[] buffer;
	} else {
		std::cout << "file could not be read" << std::endl;
	}
	if(points.size() == 0) {
		std::cout << "no points added!" << std::endl;
	}

	return points;
}
/*
Sked read( std::istream & is )
{
    // populate tree structure pt
    using boost::property_tree::ptree;
    ptree pt;
    read_xml(is, pt);

    // traverse pt
    Sked ans;
    BOOST_FOREACH( ptree::value_type const& v, pt.get_child("sked") ) {
        if( v.first == "flight" ) {
            Flight f;
            f.carrier = v.second.get<std::string>("carrier");
            f.number = v.second.get<unsigned>("number");
            f.date = v.second.get<Date>("date");
            f.cancelled = v.second.get("<xmlattr>.cancelled", false);
            ans.push_back(f);
        }
    }

    return ans;
}
 */