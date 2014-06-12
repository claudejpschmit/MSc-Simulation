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
    ct = (int)parameters[3];
    ct1 = (int)parameters[3];
    seed = (int)parameters[4];
    annealing = (int)parameters[5];
    annealMult = parameters[6];
    tMax = parameters[7];
    e0 = parameters[8];
    kb = parameters[9];

    // read interactionfile
    // initialize matrix
    interaction.resize(nPolymerase, vector<double>(nPolymerase, 0));
    //read data from file
    for (int i = 0; i < nPolymerase; ++i)
        for (int j = 0; j < nPolymerase; ++j)
            interactionFile >> interaction[i][j];
} 

void Parser::setParams(int* nsteps, int* max_index, int* nPolymerase,
            int* seed, int* ct, int* ct1, double* tMax, double* annealMult,
            double* e0, double* kb, bool* annealing)
{
    *nsteps = this->nsteps;
    *max_index = this->max_index; 
    *nPolymerase = this->nPolymerase;
    *seed = this->seed;
    *ct = this->ct;
    *ct1 = this->ct1;
    *tMax = this->tMax;
    *annealMult = this->annealMult;
    *e0 = this->e0;
    *kb = this->kb;
    *annealing = this->annealing;

}

void Parser::setInteraction(vector<vector<double>> *interaction)
{
    *interaction = this->interaction;
}
