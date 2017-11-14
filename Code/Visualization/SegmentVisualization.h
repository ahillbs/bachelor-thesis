//
// Created by gharg on 23.06.17.
//

#ifndef BA_SEGMENTVISUALIZATION_H
#define BA_SEGMENTVISUALIZATION_H


#include <utility>
#include <vector>
#include <map>
#include <set>
#include "../Capturing/Header/Segment.h"
#include "simple_svg_1.0.0.hpp"

using namespace svg;

class SegmentVisualization {
private:
    std::vector< std::vector<Vector2D>> segmentToVectorMap;
    std::set<std::pair<Vector2D,Vector2D>> edges;
    std::vector<Vector2D> takenVertices;
    Document doc;
    int xBound,yBound;
public:
    SegmentVisualization(const char* name, int xBound, int yBound);
    void setSegments(std::vector< std::vector<Vector2D>> &segmentToVectorMap);
    void setTakenEdges(std::set<std::pair<Vector2D,Vector2D>> &edges);
    void setTakenVertices(std::vector<Vector2D> &vertices);
    void drawSvg();
    void saveDocument();

};


#endif //BA_SEGMENTVISUALIZATION_H
