//
//  Point.hpp
//  Section3
//
//  Created by Benjamin Jaderberg on 17/01/2017.
//  Copyright © 2017 BenJad. All rights reserved.
//

#ifndef Point_hpp
#define Point_hpp

class Point{
    
    double xval, yval;
    
public:
    Point(double x=0.0, double y=0.0);
    
    double x();
    double y();
    void setX(double x);
    void setY(double y);
    void setPoint(double x, double y);
    double dist(Point);
    double dot(Point);
    double length();
    Point add(Point);
    Point sub(Point);
    void move(double, double);
    void copy(Point);
    void print();
    
    
    
};


#endif /* Point_hpp */
