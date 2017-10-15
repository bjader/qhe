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
void writeMCToFile(vector<vector<double>> list, string name);
void writeToFile(vector<Point> list, string name);

//Task 3 methods
complex<double> createWaveFunction(vector<Point> points);
complex<double> createWaveFunctionTerm(vector<Point> points, int option);
vector<double> calcReducedPsi(vector<Point> points);

#endif /* methods_hpp */


