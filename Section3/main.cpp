//
//  main.cpp
//  Section3
//
//  Created by Benjamin Jaderberg on 17/01/2017.
//  Copyright © 2017 BenJad. All rights reserved.
//

#include <iostream>

#include "Point.hpp"
#include "methods.hpp"

#include <cmath>
#include <complex>
#include <random>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//Global constants
const complex<double> i(0.0,1.0);
double L = 1;

//Global objects
random_device rd;
//Mersenne Twister algorithm, default seed
mt19937 gen(rd());
uniform_real_distribution<> dis(0,1); // Define our distribution as between 0 and 1


int main(int argc, const char * argv[]) {
    
    vector<Point> points;
    
    double N = 5;
    
    for (int i=0; i<N; i++) {
        double rand_x = (dis(gen) * L) - (L/2);
        double rand_y = (dis(gen) * L) - (L/2);
        
        points.push_back(Point(rand_x,rand_y));
    }
    
    complex<double> Psi = createWaveFunction(points);
  
}
