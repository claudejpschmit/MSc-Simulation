#pragma once

#include <iostream>
#include <vector>

using namespace std;

double dist(double r1[3], double r2[3]); 
void streamfile(int walkSize, vector<vector<double>> z, vector<double> polymerase);
void metropolis(double ti, int seed, vector<vector<double>> zn);
