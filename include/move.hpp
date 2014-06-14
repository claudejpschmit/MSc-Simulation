#pragma once

#include <iostream>
#include <vector>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

using namespace std;

class Move {

public: 
    Move(int seed, int walkSize);
    void pivot(vector<vector<double>>* zn, vector<vector<double>> z);
    void krank(vector<vector<double>>* zn, vector<vector<double>> z);

private:
    boost::mt19937 rng;
    boost::random::uniform_01<> rand01;
    int walkSize;
    void dir(double costh, double sinth, double cospsi, double sinpsi, 
            double cosom, double sinom, double rot[3][3]);
    void inve(double costh, double sinth, double cospsi, double sinpsi, 
            double rot[3][3]);

};
