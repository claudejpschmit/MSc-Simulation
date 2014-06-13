#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "parser.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

using namespace std;

int main(int argc, char *argv[]) {

    if (argc < 3) {
        cout << "Run this program as: " << argv[0] 
            << " interaction_file parameter_file" << endl; 
        return -1;
    }

    // initializing Programm parameters. 
    int nsteps, walkSize, max_index, nPolymerase, seed, ct, loopSize;
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
    files.setParams(&nsteps, &max_index, &nPolymerase, &seed, &ct, 
            &loopSize, &tMax, &annealMult, &e0, &kb, &annealing, &walkSize);

    //initialize RNG
    boost::mt19937 rng;
    rng.seed(seed);
    boost::random::uniform_real_distribution<> uniDistr(0.0, 1.0);
    //define RW size; walkSize == L;
    // initialize RW
    vector<vector<double>> z(walkSize, vector<double>(3, 0));
    vector<vector<double>> z_new = z;

    
    //output interaction parameters. 
    //check for interaction to be read correctly
    ofstream output_ip("interactionparameters.dat");
    //output polymerase positions
    //check if polymerase is correctly initialized
    ofstream output_pp("polymerasepositions.dat");

    for (int i = 0; i < walkSize; ++i) {
        if (polymerase[i] != 0) {
            output_pp << i << endl;
        }
        for (int j = 0; j < walkSize; ++j) {
            if (polymerase[i] + polymerase[j] == 2) {
                output_ip << i << " " << j << " " << interaction[i][j] << endl;
            }
        }
    }
    output_ip.close();
    output_pp.close();

    // defines energy output file
    ofstream output_en("energy.dat");

    // ------- Index Loop ------- //
    for (int index = 0; index < max_index; ++index) {
        double t = tMax;
        for (int i = 0; i < walkSize; ++i) {
            z[i][0] = 0.0;
            z[i][1] = 0.0;
            z[i][2] = i;
        }
    
        double en = 0;
        //int nr = 0; was in the c code but doesn't do anything

        // ----- Main Loop ----- //
        for (int n = 0; n < nsteps; ++n) {
            if (annealing)
                if (n % ct == 0)
                    t *= annealMult;
        
            double ti = 1.0 / t;

            double ev = en;
            //TODO: use marsenne twister here, random between 0 and 1
            double rand01 = uniDistr(rng);
            
            if (rand01 < 0.5) 
                en += 1; //TODO: Delete this line
                //pivot(seed, z_new);
            else
                en -=1; //TODO: Delete this line
                //krank(seed, z_new);

            //metropolis(ti, seed, z_new);

            if (n % ct == 0 || n == 0) {
                //streamfile(n); // not sure what this does
                output_en << t << " " << en << endl;
            }
        }
    }
    output_en.close();
        
    return 0;
}
