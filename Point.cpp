//
//  Point.cpp
//  Section3
//
//  Created by Benjamin Jaderberg on 17/01/2017.
//  Copyright © 2017 BenJad. All rights reserved.
//

#include "Point.hpp"
#include <iostream>
#include <math.h>
#include <stdio.h>

//Constructor
Point::Point(double x, double y) {
    xval = x;
    yval = y;
}

//Getter method
double Point::x() { return xval; }
double Point::y() { return yval; }

//Setter method
void Point::setX(double x) { xval = x; }
void Point::setY(double y) { yval = y; }
void Point::setPoint(double x, double y) { xval = x; yval = y;}


// Distance to another point.  Pythagorean thm.
double Point::dist(Point other) {
    double xd = xval - other.xval;
    double yd = yval - other.yval;
    return sqrt(xd*xd + yd*yd);
}

//Scalar product between two points
double Point::dot(Point other) {
    
    double dot_product = (xval*other.xval) + (yval*other.yval);
    
    return dot_product;
}

//Length of a vector point
double Point::length() {
    return sqrt(xval*xval + yval*yval);
}

// Add or subtract two points.
Point Point::add(Point b)
{
    return Point(xval + b.xval, yval + b.yval);
}
Point Point::sub(Point b)
{
    return Point(xval - b.xval, yval - b.yval);
}

// Move the existing point.
void Point::move(double a, double b)
{
    xval += a;
    yval += b;
}

//Make a copy of another point
void Point::copy(Point other) {
    xval = other.xval;
    yval = other.yval;
}

// Print the point on the stream.  The class ostream is a base class for output streams of various types.
void Point::print()
{
    std::cout << "(" << xval << "," << yval << ")";
}


