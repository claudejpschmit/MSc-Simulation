#include "helper.hpp"

#include <cmath>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>


double dist(double r1[3], double r2[3]) {
double distance = pow(r1[0] - r2[0],2.0) + pow(r1[1] - r2[1],2.0) + 
    pow(r1[2] - r2[2],2.0);

return sqrt(distance);
}

void streamfile(int step)
{

}
void metropolis(double ti, int seed, vector<vector<double>> zn)
{

}
