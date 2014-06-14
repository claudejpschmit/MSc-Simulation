#include "helper.hpp"

#include <cmath>
#include <fstream>

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
