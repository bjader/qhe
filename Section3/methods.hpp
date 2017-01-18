//
//  methods.hpp
//  Section3
//
//  Created by Benjamin Jaderberg on 17/01/2017.
//  Copyright © 2017 BenJad. All rights reserved.
//

#include "Point.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <complex>

using namespace std;

#ifndef methods_hpp
#define methods_hpp

//Task 2 methods
void writeMCToFile(vector<Point> list, vector<double> stdev, string name);
double calcProbability (complex<double> psi1, complex<double> psi2);

//Task 3 methods
complex<double> createWaveFunction(vector<Point> points);


#endif /* methods_hpp */


