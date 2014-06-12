#pragma once

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class Parser {

public: 
    Parser(ifstream &interactionFile, ifstream &parameterFile);
    void setParams(int* nsteps, int* max_index, int* nPolymerase,
            int* seed, int* ct, int* ct1, double* tMax, double* annealMult,
            double* e0, double* kb, bool* annealing);
    void setInteraction(vector<vector<double>>* interaction);
private: 
    int nsteps, max_index, nPolymerase, seed, ct, ct1;
    double tMax, annealMult, e0, kb;
    bool annealing;
    vector<vector<double>> interaction;


};
