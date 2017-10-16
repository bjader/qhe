//
//  methods.cpp
//  Section3
//
//  Created by Benjamin Jaderberg on 17/01/2017.
//  Copyright © 2017 BenJad. All rights reserved.
//

#include "methods.hpp"

//Method to write MC generated data to file, a vector of Points and a vector of stdev doubles
void writeMCToFile(vector<vector<double>> list, string name) {
    ofstream myfile("/Users/benjaminjaderberg/git/msci-project/qhe/data/" + name);
    if (myfile.is_open()) {
        for (int i=0; i<list.size(); i++) {
            myfile << (list[i])[0] << "," << (list[i])[1] << "," << (list[i])[2] << "\n";
        }
        myfile.close();
    }
    else { cout << "Writing to file failed";}
}

//Method to write a vector of points to comma seperated file
void writeToFile(vector<Point> list, string name) {
    ofstream myfile("/Users/benjaminjaderberg/Desktop/4th_Year/MSci_Project/Section3/" + name);
    if (myfile.is_open()) {
        for (Point p : list) {
            myfile << p.x() << "," << p.y() << "\n";
        }
        myfile.close();
    }
    else { cout << "Writing to file failed";}
}

//Method to generate both ln(|Psi|) and the complex component of Psi
//Used instead of createWaveFunction() when Psi > 10^300 and becomes classified as 'inf'


vector<double> calcReducedPsi(vector<Point> points, int m) {
    
    double logPsi = 0;
    double Phi = 0;
    vector<complex<double>> zPoints;
    
    for (Point p : points) {
        zPoints.push_back(complex<double>(p.x(),p.y()));
    }
        
    //For every point in the system
    for (int i=0; i<zPoints.size(); i++) {
        
        complex<double> z_i = zPoints[i];
        logPsi -= norm(z_i);
        
        for (int j=i+1; j<zPoints.size(); j++) {
            
            complex<double> z_j = zPoints[j];
            logPsi += m*log(abs(z_i-z_j));
            
            Phi += m*arg(z_i-z_j);
            
        }
    }
    
    vector<double> output = {logPsi, Phi};
    
    return output;
    
}

