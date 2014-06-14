#include "metropolis.hpp"

#include <cmath>

Metropolis::Metropolis(int walkSize, double e0, double kb, 
        vector<vector<double>> interaction, vector<double> polymerase)
    :
        walkSize(walkSize),
        kb(kb),
        e0(e0)
{
    this->interaction = interaction;
    this->polymerase = polymerase;
}

void Metropolis::step(double* en, double ev, double ti, vector<vector<double>>* z, 
            vector<vector<double>> zn, double rand01)
{
    double r0 = 0.95;
    double r1 = 1.4;
    double rag1;
    double v1[3];
    double v2[3];

    vector<vector<double>> tang(walkSize, vector<double> (3, 0));
    for (int i = 0; i < walkSize - 1; ++i)
        for (int j = 0; j < 3; ++ j)
            tang[i][j] = zn[i + 1][j] - zn[i][j];
    double energy = 0.0;

    for (int i = 0; i < walkSize - 2; ++i)
        energy -= kb * (tang[i][0] * tang[i + 1][0]
                + tang[i][1] * tang[i + 1][1]
                + tang[i][2] * tang[i + 1][2]);

    for (int i = 0; i < walkSize - 1; ++ i) {
        for (int j = i + 1; j < walkSize; ++j) {
            if (j != i + 1) {
                if (polymerase[i] + polymerase[j] == 2) {
                    for (int k = 0; k < 3; ++k) {
                        v1[k] = zn[j][k];
                        v2[k] = zn[i][k];
                    }
                    rag1 = dist(v1, v2);
                    if (rag1 < r0) {
                        energy -= 100000.0;
                    }
                    if (rag1 >= r0 && rag1 <= r1) {
                        if (interaction[i][j] == 1) {
                            energy += e0;
                        }
                    }
                }
            }
        }
    }

    double  asc = ti * (energy - ev);

    if (log(rand01) < asc) 
        *z = zn;

    if (log(rand01) >= asc)
        energy = ev;

    *en = energy;

}

double Metropolis::dist(double r1[3], double r2[3]) 
{
    double distance = pow(r1[0] - r2[0],2.0) + pow(r1[1] - r2[1],2.0) + 
        pow(r1[2] - r2[2],2.0);

    return sqrt(distance);
}
