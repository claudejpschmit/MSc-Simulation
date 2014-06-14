#include "helper.hpp"

#include <cmath>
#include <fstream>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>


double dist(double r1[3], double r2[3]) {
double distance = pow(r1[0] - r2[0],2.0) + pow(r1[1] - r2[1],2.0) + 
    pow(r1[2] - r2[2],2.0);

return sqrt(distance);
}

void streamfile(int walkSize, vector<vector<double>> z, vector<double> polymerase)
{
    ofstream out1("configuration.dat");
    ofstream out2("polymerases.dat");
    ofstream out3("distancecheck.dat");

    for (int i = 0; i < walkSize; ++i) {
        out1 << z[i][0] << " " << z[i][1] << " " << z[i][2] << endl;
        if (polymerase[i] == 1)
            out2 << z[i][0] << " " << z[i][1] << " " << z[i][2] << endl;
        if (i != walkSize - 1) { 

            double diff1 = z[i + 1][0] - z[i][0];
            diff1 *= diff1;
            double diff2 = z[i + 1][1] - z[i][1];
            diff2 *= diff2;
            double diff3 = z[i + 1][2] - z[i][2];
            diff3 *= diff3;

            out3 << i << " " << sqrt(diff1 + diff2 + diff3) << endl;
        }
    }

    out1.close();
    out2.close();
    out3.close();
}
void metropolis(double ti, int seed, vector<vector<double>> zn)
{

}
