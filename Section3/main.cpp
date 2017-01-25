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
//const complex<double> i(0.0,1.0);
double L = 5;

//Global objects
random_device rd;
//Mersenne Twister algorithm, default seed
mt19937 gen(rd());
uniform_real_distribution<> dis(0,1); // Define our distribution as between 0 and 1

//Method to "burn in" a random system to more accurately represent the true distribution
//Runs same Metropolis algorithm as runMetropolis but without any calculations or swaps
void runBurnIn (vector<Point> &R1, vector<Point> &R2, double dr, int num_iterations) {
    
    random_device rd;
    //Mersenne Twister algorithm, default seed
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0,1);
    
    for (int i = 0; i<num_iterations; i++) {
        
        //Calculate probability of current system
        complex<double> psi1 = createWaveFunction(R1);
        complex<double> psi2 = createWaveFunction(R2);
        double p = calcJointProb(psi1, psi2);
        
        //Create empty test system
        vector<Point> R3;
        vector<Point> R4;
        
        //Populate new system R3 with original particles from R1, plus a random displacement of length dr
        for (Point p1 : R1) {
            double phi = (2*M_PI) * dis(gen);
            double x = p1.x() + (dr * cos(phi));
            double y = p1.y() + (dr * sin(phi));
            
            //If a particle moves outside the box, instead make it reappear out the opposite side
            if (x > (L/2)) { x = x - L;}
            if (x < (-L/2)) {x = x + L;}
            if (y > (L/2)) { y = y - L;}
            if (y < (-L/2)) { y = y + L;}
            
            R3.push_back(Point(x,y));
        }
        
        //Repeat for R4/R2
        for (Point p2 : R2) {
            double phi = (2*M_PI) * dis(gen);
            double x = p2.x() + (dr * cos(phi));
            double y = p2.y() + (dr * sin(phi));
            
            if (x > (L/2)) { x = x - L;}
            if (x < (-L/2)) {x = x + L;}
            if (y > (L/2)) { y = y - L;}
            if (y < (-L/2)) { y = y + L;}
            
            R4.push_back(Point(x,y));
            
        }
        
        //Calculate probability of new system
        complex<double> psi3 = createWaveFunction(R3);
        complex<double> psi4 = createWaveFunction(R4);
        double p_new = calcJointProb(psi3,psi4);
        
        //Accept in accordance to Hastings-Metropolis method
        double lambda = min(p_new/p,1.0);
        double alpha = dis(gen);
        
        //If accept the new configuration
        if (alpha < lambda) {
            
            R1 = R3;
            psi1 = psi3;
            R2 = R4;
            psi2 = psi4;
        }
        
    }
}


//Method using a Metropolis algorithm to iterate the two systems over random moves
void runMetropolis (vector<Point> rPoints1, vector<Point> rPoints2, double dr, int num_iterations) {
    
    //Define some counters
    int accepted = 0;
    int rejected = 0;
    int swapped = 0;
    
    //Set up array to keep track of integral summations for a set number of Z values
    vector<complex<double>> sum_gx;
    vector<complex<double>> sum_g2x;
    for (int a=0; a<1; a++) {
        sum_gx.push_back(0);
        sum_g2x.push_back(0);
    }
    
    //Define current system
    vector<Point> R1(rPoints1);
    vector<Point> R2(rPoints2);
    
    for (int i = 0; i<num_iterations; i++) {
        
        //Calculate probability of current system
        complex<double> psi1 = createWaveFunction(rPoints1);
        complex<double> psi2 = createWaveFunction(rPoints2);
        double p = calcJointProb(psi1, psi2);
        
        //Create empty test system
        vector<Point> R3;
        vector<Point> R4;
        
        //Populate new system R3 with original particles from R1, plus a random displacement of length dr
        for (Point p1 : R1) {
            double phi = (2*M_PI) * dis(gen);
            double x = p1.x() + (dr * cos(phi));
            double y = p1.y() + (dr * sin(phi));
            
            //If a particle moves outside the box, instead make it reappear out the opposite side
            if (x > (L/2)) { x = x - L;}
            if (x < (-L/2)) {x = x + L;}
            if (y > (L/2)) { y = y - L;}
            if (y < (-L/2)) { y = y + L;}
            
            R3.push_back(Point(x,y));
        }
        
        //Repeat for R4/R2
        for (Point p2 : R2) {
            double phi = (2*M_PI) * dis(gen);
            double x = p2.x() + (dr * cos(phi));
            double y = p2.y() + (dr * sin(phi));
            
            if (x > (L/2)) { x = x - L;}
            if (x < (-L/2)) {x = x + L;}
            if (y > (L/2)) { y = y - L;}
            if (y < (-L/2)) { y = y + L;}
            
            R4.push_back(Point(x,y));
            
        }
        
        //Calculate probability of new system
        complex<double> psi3 = createWaveFunction(R3);
        complex<double> psi4 = createWaveFunction(R4);
        double p_new = calcJointProb(psi3,psi4);
        
        //Accept in accordance to Hastings-Metropolis method
        double lambda = min(p_new/p,1.0);
        double alpha = dis(gen);
        
        //If accept the new configuration
        if (alpha < lambda) {
            
            accepted += 1;
            complex<double> g;
            vector<int> iListR3;
            vector<int> iListR4;
            
            //The new moved system becomes our current system
            R1 = R3;
            psi1 = psi3;
            R2 = R4;
            psi2 = psi4;
            
            //We will keep using R3 and R4 as the system for our swap calculation
                
            //Find the ID of particles within subsystem A (ie. particles left of the y axis)
            for (int j=0; j<R3.size(); j++) {
                if (R3[j].x() < 0) {
                    iListR3.push_back(j);
                }
                if (R4[j].x() < 0) {
                    iListR4.push_back(j);
                    
                }
            }
            //If the number of particles in R3_A = R4_A we can perform the swap
            if (iListR3.size() == iListR4.size()) {
                swapped += 1;
                for (int d=0; d<iListR3.size(); d++) {
                    
                    //Find the index integer of the particles we are swapping
                    int id_1 = iListR3[d];
                    int id_2 = iListR4[d];
                        
                    //Swap the particles by copying each others coordinates
                    //At this point R2 is still an exact copy of R4. So in fact we copy R2_A into R3_A to avoid creating temporary copy variables
                    //And we do the same copying R1_A into R4_A
                    R3[id_1].copy(R2[id_2]);
                    R4[id_2].copy(R1[id_1]);
                }
                //Recalculate the wavefunctions of R3 and R4 (after the swap)
                    
                psi3 = createWaveFunction(R3);
                psi4 = createWaveFunction(R4);
                g = (psi3 * psi4)/(psi1 * psi2);
                    
            }
            //If we didnt perform the swap, dont add to the integral (dirac delta factor)
            else {
                g = 0;
            }
            sum_gx[0] += g;
            sum_g2x[0] += pow(g,2);
                
            iListR3.clear();
            iListR4.clear();
            
        }
        //Else reject the new configuration
        else {
            rejected +=1;
        }
        
    }
    //Calculate the error of each S2 point
    vector<double> normIntegral;
    for (int b=0; b<sum_gx.size(); b++) {
        normIntegral.push_back((sum_gx[b].real() / accepted));
    }
    vector<double> normSquaredIntegral;
    for (int i=0; i<sum_g2x.size(); i++) {
        normSquaredIntegral.push_back(sum_g2x[i].real() / accepted);
    }
    
    vector<double> stdevS2;
    for (int i=0; i<sum_gx.size(); i++) {
        stdevS2.push_back((pow(normSquaredIntegral[i] - pow(normIntegral[i],2),0.5))/(pow(accepted,0.5)*normIntegral[i]));
    }
    vector<double> s2;
    
    //Print properties of run
    cout << endl << "Number of accepted moves: " << accepted;
    cout << endl << "Number of rejected moves: " << rejected;
    cout << endl << "Acceptance rate: " << (double(accepted)/(accepted + rejected))*100 << "%";
    cout << endl << "Number of swaps: " << swapped;
    
    //Print normalised integral and calculate S2 for each one
    cout << endl << "Normalised integrals: ";
    for (double d : normIntegral) {
        cout << d << " ";
        s2.push_back(-log(d));
    }
    
    //Print S2 with the standard deviation
    cout << endl << "S2 values: ";
    for (int i = 0; i<s2.size(); i++) {
        cout << s2[i] << "+/-" << stdevS2[i] << endl;
    }
}


int main(int argc, const char * argv[]) {
    
    vector<Point> R1;
    vector<Point> R2;
    
    double N = 5;
    
    for (int i=0; i<N; i++) {
        
        //Populate R1 with particles randomly placed within radius of 3-sigma from origin
        double sigma = 1;
        
        double phi1 = (2*M_PI) * dis(gen); //Random angle in polar coordinates
        double radius1 = (dis(gen) * 3 * sigma); //Random radius between 0 -> 3-sigma
        
        double x1 =  radius1 * cos(phi1);
        double y1 = radius1 * sin(phi1);
        R1.push_back(Point(x1,y1));
        
        //Repeat for R2
        
        double phi2 = (2*M_PI) * dis(gen); //Random angle in polar coordinates
        double radius2 = (dis(gen) * 3 * sigma); //Random radius between 0 -> 3-sigma
        
        double x2 =  radius2 * cos(phi2);
        double y2 = radius2 * sin(phi2);
        R2.push_back(Point(x2,y2));
    }
    
    runBurnIn(R1, R2, 0.1, 100000);

    runMetropolis(R1, R2, 0.005, 100000);
    
}
