#include "move.hpp"

#include <cmath>

#define PI 3.14159265359

Move::Move(int seed, int walkSize) 
    :
        walkSize(walkSize)
{
    rng.seed(seed);    
}

void Move::pivot(vector<vector<double>>* zn, vector<vector<double>> z) 
{
    vector<vector<double>> buff = *zn;
    double ver = 1;
    if (rand01(rng) < 0.5) ver = -1;
    double beta = rand01(rng) * PI * ver;   //max angle pi

    double teta = acos(rand01(rng)); 

    ver = 1;
    if (rand01(rng) < 0.5) ver = -1;
    double fi = rand01(rng) * PI * ver;   //max angle pi

    int nm = floor(rand01(rng) * (walkSize - 1));

    int maxI = walkSize - nm;
    vector<vector<double>> pos(walkSize, vector<double>(3, 0));
    vector<vector<double>> pos1(walkSize, vector<double>(3, 0));

    for (int i = 0; i < maxI; ++i)
        for (int j = 0; j < 3; ++j)
            pos[i][j] = z[nm + i][j] - z[nm][j];

    double rot[3][3];
    dir(cos(teta), sin(teta), cos(fi), sin(fi), cos(beta), sin(beta), rot);

    for (int i = 0; i < maxI; ++i) {
        for (int j = 0; j < 3; ++j) {
            double sum = 0;
            for (int k = 0; k < 3; ++k) {
                sum += rot[j][k] * pos[i][k];
            }
            pos1[i][j] = sum;
        }
    }

    inve(cos(teta), sin(teta), cos(fi), sin(fi), rot);

    for (int i = 0; i < maxI; ++i) {
        for (int j = 0; j < 3; ++j) {
            double sum = 0;
            for (int k = 0; k < 3; ++k) {
                sum += rot[j][k] * pos1[i][k];
            }
            pos[i][j] = sum;
        }
    }

    for (int i = 0; i < maxI; ++i)
        for (int j = 0; j < 3; ++j)
            buff[nm + i][j] = pos[i][j] + z[nm][j];

    if (nm != 0)
        for (int i = 0; i < nm; ++i)
            for (int j = 0; j < 3; ++j)
                buff[i][j] = z[i][j];
    
    *zn = buff;
}

void Move::krank (vector<vector<double>>* zn, vector<vector<double>> z) 
{
    vector<vector<double>> pos(walkSize, vector<double>(3, 0));
    vector<vector<double>> pos1(walkSize, vector<double>(3, 0));
    vector<vector<double>> buff = *zn;

Choose: int nm = floor(rand01(rng) * walkSize);
    if (nm >= walkSize - 2 || nm == 0) goto Choose;

    int nm1 = floor(rand01(rng) * 100 + 2);
    if ((nm + nm1) >= walkSize - 1) nm1 = walkSize - 1 - nm;

    int nf = nm + nm1;
    int ver = 1;
    if (rand01(rng) < 0.5) ver = - 1;
    double beta = rand01(rng) * PI * ver;           // max angle pi
    /*
    double base[3];
    for (int i = 0; i < 3; ++i)
        base[i] = z[nf][i];
*/

    double c[3];
    for (int i = 0;i < 3; ++i)
        c[i] = z[nm][i];

    for (int i = 0; i < nm1 + 1; ++i)
        for (int j = 0; j < 3; ++j)
            pos[i][j] = z[nm + i][j] - c[j];

    double rag = pow(pos[nm1][0], 2.0) + pow(pos[nm1][1], 2.0) + 
        pow(pos[nm1][2], 2.0);
    rag = sqrt(rag);
    double costh = pos[nm1][2] / rag;
    double sinth = sqrt(1.0 - costh * costh);
    double cospsi, sinpsi, cosom, sinom;

    if (sinth < 0.0001) goto ReAssembleChain;
    cospsi = pos[nm1][0] / (rag * sinth);
    sinpsi = pos[nm1][1] / (rag * sinth);
    cosom = cos(beta);
    sinom = sin(beta);
    double rot[3][3];
    dir(costh, sinth, cospsi, sinpsi, cosom, sinom, rot);  

    for (int i = 0; i < nm1 + 1; ++i){
        for (int j = 0; j < 3; ++j){
            double sum = 0.0;
            for (int k = 0; k < 3; ++k){
                sum += rot[j][k] * pos[i][k];
            }
            pos1[i][j] = sum;
        }
    }

    inve(costh,sinth,cospsi,sinpsi,rot);

    for (int i = 0; i < nm1 + 1; ++i){
        for (int j = 0; j < 3; ++j){
            double sum = 0.0;
            for (int k = 0; k < 3; ++k){
                sum += rot[j][k] * pos1[i][k];
            }
            pos[i][j] = sum;
        }
    }

ReAssembleChain: for(int i = 0; i < nm1 + 1; ++i){
                        for(int j = 0; j < 3; ++j){
                            buff[nm + i][j] = pos[i][j] + c[j];
                        }  
                    }

    if (nm != 0) 
        for (int i = 0; i < nm; ++i)
            for (int j = 0; j < 3; ++j)
                buff[i][j] = z[i][j];

    if (nf != walkSize - 1)
        for(int i = nf + 1; i <= walkSize - 1; ++i)
            for(int j = 0; j < 3; ++j)
                buff[i][j] = z[i][j];

    /*
    for (int i = 0; i < 3; ++i){
        if (abs(base[i] - buff[nf][i]) > 0.01){
            fprintf(stderr, "aarrghhh!!!!\n");
            fprintf(stderr, "%d %lf\n", i, base[i] - buff[nf][i]);
        }
    }
    */

    *zn = buff;
}  

void Move::dir(double costh, double sinth, double cospsi, double sinpsi, 
        double cosom, double sinom, double rot[3][3])
{
    rot[0][0] = cospsi * costh * cosom + sinom * sinpsi;
    rot[0][1] = cosom * costh * sinpsi - sinom * cospsi;
    rot[0][2] = - cosom * sinth;
    rot[1][0] = sinom * cospsi * costh - sinpsi * cosom;
    rot[1][1] = sinpsi * sinom * costh + cosom * cospsi;
    rot[1][2] = - sinom * sinth;
    rot[2][0] = cospsi * sinth;
    rot[2][1] = sinpsi * sinth;
    rot[2][2] = costh;
}
void Move::inve(double costh, double sinth, double cospsi, double sinpsi, 
        double rot[3][3])
{
    rot[0][0] = cospsi * costh;
    rot[0][1] = - sinpsi;
    rot[0][2] = cospsi * sinth;
    rot[1][0] = sinpsi * costh;
    rot[1][1] = cospsi;
    rot[1][2] = sinpsi * sinth;
    rot[2][0] = - sinth;
    rot[2][1] = 0;
    rot[2][2] = costh;
}
