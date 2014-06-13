#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "parser.hpp"

using namespace std;

int main(int argc, char *argv[]) {

    if (argc < 3) {
        cout << "Run this program as: " << argv[0] 
            << " interaction_file parameter_file" << endl; 
        return -1;
    }

    // initializing Programm parameters. 
    int nsteps, walkSize, max_index, nPolymerase, seed, ct, ct1, loopSize;
    double tMax, annealMult, e0, kb;
    bool annealing;
    vector<vector<double>> interaction;
    vector<double> polymerase;

    ifstream parameterFile(argv[2]);
    ifstream interactionFile(argv[1]);
    Parser files(interactionFile, parameterFile);
    parameterFile.close();
    interactionFile.close();
    files.setInteraction(&polymerase, &interaction);
    files.setParams(&nsteps, &max_index, &nPolymerase, &seed, &ct, &ct1, 
            &loopSize, &tMax, &annealMult, &e0, &kb, &annealing, &walkSize);


    // Checks! 
    // parameters initialize correctly
    /*
    cout << nsteps << endl;
    cout << max_index << endl;
    cout << nPolymerase << endl;
    cout << seed << endl;
    cout << ct << endl;
    cout << ct1 << endl;
    cout << loopSize << endl;
    cout << tMax << endl;
    cout << annealMult << endl;
    cout << e0 << endl;
    cout << kb << endl;
    cout << annealing << endl;
    */    

    // interaction initializes correctly
   /* 
    for (int i = 0; i < nPolymerase; ++i) {
        for (int j = 0; j < nPolymerase; ++j) {
            cout << interaction[i][j] << "  ";
        }
        cout << endl;
    }
    */

    //define RW size; walkSize == L;
    // initialize RW
    vector<vector<double>> z(walkSize, vector<double>(3, 0));
    vector<vector<double>> z_new = z;

    
    //output interaction parameters. 
    ofstream output_ip("interactionparameters.dat");
    for (int i = 0; i < walkSize; ++i) {
        for (int j = 0; j < walkSize; ++j) {
            if (polymerase[i] + polymerase[j] == 2) {
                output_ip << i << " " << j << " " << interaction[i][j] << endl;
            }
        }
    }
    output_ip.close();

    
    return 0;
}
