#pragma once

#include <iostream>
#include <vector>

using namespace std;

class Metropolis {

public:
    Metropolis(int walkSize, double e0, double kb,
            vector<vector<double>> interaction, vector<double> polymerase);
    void step(double* en, double ev, double ti, vector<vector<double>>* z, 
            vector<vector<double>> zn, double rand01);
private:
    double dist(double r1[3], double r2[3]);
    int walkSize;
    double kb, e0;
    vector<vector<double>> interaction;
    vector<double> polymerase;
};
