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
    int nsteps, max_index, nPolymerase, seed, ct, ct1;
    double tMax, annealMult, e0, kb;
    bool annealing;
    vector<vector<double>> interaction;

    ifstream parameterFile(argv[2]);
    ifstream interactionFile(argv[1]);
    Parser files(interactionFile, parameterFile);
    parameterFile.close();
    interactionFile.close();
    files.setInteraction(&interaction);
    files.setParams(&nsteps, &max_index, &nPolymerase, &seed, &ct, &ct1, 
            &tMax, &annealMult, &e0, &kb, &annealing);

    // Checks! 
    // parameters initialize correctly
    /*
    cout << nsteps << endl;
    cout << max_index << endl;
    cout << nPolymerase << endl;
    cout << seed << endl;
    cout << ct << endl;
    cout << ct1 << endl;
    cout << tMax << endl;
    cout << annealMult << endl;
    cout << e0 << endl;
    cout << kb << endl;
    cout << annealing << endl;
    */

    // interaction initializes correctly
    /*
    for (int j = 0; j < nPolymerase; ++j) {
        for (int i = 0; i < nPolymerase; ++i) {
            cout << interaction[i][j] << "  ";
        }
        cout << endl;
    }
    */
    return 0;
}
