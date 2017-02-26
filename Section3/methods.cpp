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
    ofstream myfile("/Users/benjaminjaderberg/Desktop/4th_Year/MSci_Project/Section3/" + name);
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


//Method to calculate the probability over two given wave function: p = |psi1|^2 * |psi2|^2
double calcJointProb (complex<double> psi1, complex<double> psi2) {
    
    return norm(psi1) * norm(psi2);
}

//Method to create a Laughlin wavefunction based off input coordinates in 2-d space
complex<double> createWaveFunction(vector<Point> points) {
    
    double scale_factor = 2;
    
    complex<double> Psi;
    vector<complex<double>> zPoints;
    
    for (Point p : points) {
        zPoints.push_back(complex<double>(p.x(),p.y()));
    }
    
    //Initialise summation variables
    complex<double> product = complex<double>(1.0,0.0);
    double sum_square = 0;
    
    //For every point in the system
    for (int i=0; i<zPoints.size(); i++) {
        
        complex<double> z_i = zPoints[i];
        
        for (int j=0; j<zPoints.size(); j++) {
            
            complex<double> z_j = zPoints[j];
            complex<double> product_term = (z_i - z_j);
            //cout << endl << product_term;
            
            if(product_term == complex<double>(0.0,0.0)) {
            }
            else {
                product = product * (product_term/scale_factor);
            }
        
        }
        
        //Add |Z|^2 to the exponent sum
        sum_square += norm(z_i);
    }
    
    Psi = product * exp(-sum_square);

    return Psi;
}

//Method to calculate only one term of the Laughlin wave function
complex<double> createWaveFunctionTerm(vector<Point> points, int option) {
    
    double scale_factor = 2;
    
    complex<double> Psi;
    vector<complex<double>> zPoints;
    
    for (Point p : points) {
        zPoints.push_back(complex<double>(p.x(),p.y()));
    }
    
    //Initialise summation variables
    complex<double> product = complex<double>(1.0,0.0);
    double sum_square = 0;
    
    //For every point in the system
    for (int i=0; i<zPoints.size(); i++) {
        
        complex<double> z_i = zPoints[i];
        
        for (int j=0; j<zPoints.size(); j++) {
            
            complex<double> z_j = zPoints[j];
            complex<double> product_term = (z_i - z_j);
            
            if(product_term == complex<double>(0.0,0.0)) {
            }
            else {
                product = product * (product_term/scale_factor);
            }
            
        }
        
        //Add |Z|^2 to the exponent sum
        sum_square += norm(z_i);
    }
    
    complex<double> exp_term = exp(-sum_square);
    
    if (option) {
        return exp_term;
    }
    else {
        return product;
    }
}

