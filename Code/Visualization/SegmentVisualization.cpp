//
// Created by gharg on 23.06.17.
//

#include "SegmentVisualization.h"

SegmentVisualization::SegmentVisualization(const char *name, int xBound, int yBound) : doc(name,svg::Dimensions(xBound+5,yBound+5)), xBound(xBound), yBound(yBound) {

}

void SegmentVisualization::setTakenEdges(std::set<std::pair<Vector2D, Vector2D>> &edges) {
    this->edges = edges;
}

void SegmentVisualization::setSegments(std::vector< std::vector<Vector2D>> &segmentToVectorMap) {
    this->segmentToVectorMap = segmentToVectorMap;
}


void SegmentVisualization::setTakenVertices(std::vector<Vector2D> &vertices) {
    this->takenVertices = vertices;
}


void SegmentVisualization::drawSvg() {
    //LineChart chart(Dimensions(5,5));
    Point offset(5,5);
    double circleDiameter = 3;
    std::vector<Polyline> polyLines;
    std::set<Circle> circles;
    std::vector<Line> lines;
    for( std::vector<Vector2D> secVecs : segmentToVectorMap) {
        //Polyline line(Stroke(1, Color::Blue));
		Point startPoint(CGAL::to_double(secVecs[0].x),CGAL::to_double(secVecs[0].y));
		Point endPoint(CGAL::to_double(secVecs[secVecs.size()-1].x),CGAL::to_double(secVecs[secVecs.size()-1].y));
		Line singleLine(startPoint,endPoint,Stroke(1, Color::Blue));
        //for(Vector2D vector : secVecs) {
        //    line << Point(vector.x,vector.y);
        //}
        singleLine.offset(offset);
        //polyLines.push_back(line);
		lines.push_back(singleLine);
        doc << singleLine;
		for(Vector2D vector : secVecs) {
			Circle c(Point(CGAL::to_double(vector.x),CGAL::to_double(vector.y)),circleDiameter,Color::Black);
			c.offset(offset);
			if(circles.find(c) == circles.end())
				circles.insert(c);

		}
		/*
        for (int i = 0; i < line.points.size(); ++i) {
            Circle c(line.points[i],circleDiameter,Color::Black);
            if(circles.find(c) == circles.end())
                circles.insert(c);

        }
        */
    }

    for(std::pair<Vector2D,Vector2D> vecPair : edges) {
        Polyline line(Stroke(1, Color::Red));
        line << Point(CGAL::to_double(vecPair.first.x),CGAL::to_double(vecPair.first.y)) <<
			 Point(CGAL::to_double(vecPair.second.x),CGAL::to_double(vecPair.second.y));
        line.offset(offset);
        polyLines.push_back(line);
        doc << line;
    }

    for (unsigned int i = 0; i < takenVertices.size(); ++i) {
        Circle c(Point(CGAL::to_double(takenVertices[i].x),CGAL::to_double(takenVertices[i].y)),circleDiameter,Color::Black);
        c.offset(offset);
        auto iterator = circles.find(c);
        if(iterator != circles.end()) {
            circles.erase(c);
        }
    }
	for(Circle c : circles) {
		doc << c;
	}
	circles.clear();
	for (unsigned int i = 0; i < takenVertices.size(); ++i) {
		Circle c2(Point(CGAL::to_double(takenVertices[i].x),CGAL::to_double(takenVertices[i].y)),circleDiameter,Color::Yellow);
		c2.offset(offset);
		circles.insert(c2);
	}
	for(Circle c : circles) {
		doc << c;
	}
}

void SegmentVisualization::saveDocument() {
    doc.save();
}

