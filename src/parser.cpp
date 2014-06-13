#include "parser.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

Parser::Parser(ifstream &interactionFile, ifstream &parameterFile)
{
    // read Parameterfile
    string line;
    vector<double> parameters;
    while (!parameterFile.eof()) {
        getline(parameterFile, line);
        istringstream linestream(line);
        if (line[0] != '#' && line.length() != 0) {
            double x;
            linestream >> x;
            parameters.push_back(x);
        }
    }
    
    // define parameters
    nsteps = (int)parameters[0];
    max_index = (int)parameters[1];
    nPolymerase = (int)parameters[2];
    loopSize = (int)parameters[3];
    ct = (int)parameters[4];
    seed = (int)parameters[5];
    annealing = (int)parameters[6];
    annealMult = parameters[7];
    tMax = parameters[8];
    e0 = parameters[9];
    kb = parameters[10];

    // read interactionfile
    // initialize matrix
    walkSize = loopSize * (nPolymerase - 1) + 1;
    interaction.resize(walkSize, vector<double>(walkSize, 0));
    polymerase.resize(walkSize, 0);
    
    // initialize polymerase
    for (int i = 0; i < walkSize; ++i) {
        if (i % loopSize == 0) polymerase[i] = 1;
    }
    //read data from file
    for (int i = 0; i < walkSize; ++i) {
        if (polymerase[i] == 1) {
            for (int j = 0; j < nPolymerase; ++j) {                
                interactionFile >> interaction[i][j * loopSize];
            }
        }
    }
} 

void Parser::setParams(int* nsteps, int* max_index, int* nPolymerase,
            int* seed, int* ct, int* loopSize, double* tMax, 
            double* annealMult, double* e0, double* kb, bool* annealing, 
            int* walkSize)
{
    *nsteps = this->nsteps;
    *max_index = this->max_index; 
    *nPolymerase = this->nPolymerase;
    *seed = this->seed;
    *ct = this->ct;
    *loopSize = this->loopSize;
    *tMax = this->tMax;
    *annealMult = this->annealMult;
    *e0 = this->e0;
    *kb = this->kb;
    *annealing = this->annealing;
    *walkSize = this->walkSize;
}

void Parser::setInteraction(vector<double>* polymerase,
        vector<vector<double>> *interaction)
{
    *interaction = this->interaction;
    *polymerase = this->polymerase;
}
