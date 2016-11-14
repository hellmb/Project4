/*
 argv[1] : number of spins in each direction
 argv[2] : number of Monte Carlo cycles
 argv[3] : initial temperature
 argv[4] : decides which loop we want to enter and which files we want to write to
*/

#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <fstream>
#include "functions.h"
#include "mpi.h"

using namespace std;
using namespace arma;


int main(int argc, char* argv[]){

    if (argc <= 4){
        cout << "ERROR: Bad usage! Please see the README.mp file for command line arguments.";
        exit(1);
    }

    // convert input arguments
    int n = atoi( argv[1] );
    int mc_cycles = atoi( argv[2] );
    double T_init = atof( argv[3] );

    double beta = 1.0 / T_init;
    double w[17];

    // initialize array
    for (int dE = -8; dE <= 8; dE++){ w[dE + 8] = 0; }
    // every eight element is an exponential of energy(-8, -4, 0, 4, 8)
    for (int dE = -8; dE <= 8; dE+=4){ w[dE + 8] = exp(- beta * dE); }

    // counting accepted configurations
    int accepted_configs = 0;

    mat L = zeros<mat>(n, n);

    // make each random number be random every time we call rand()
    srand (time(NULL));

    // create a random matrix
    CreateRandomMatrix(n, L);

    // create an ordered matrix
    //CreateMatrix(n, L);

    // initialize energy and magnetization
    double E = CalculateEnergy(n, L);
    double M = CalculateMagneticMoment(L);

    // values to add to after each Monte Carlo cycle
    double energy_sum = E;
    double mag_sum = M;
    double energy_sqr = E * E;
    double mag_sqr = M * M;


    // avoid creating empty files when not running specific simulation
    if (atoi(argv[4]) == 1){

        // files to write to
        ofstream plotfile, acceptfile;

        string out1 ("../files_python/mc_cycles");
        out1 += argv[1];
        out1 += ".txt";
        plotfile.open(out1);

        string out2 ("../files_python/accepted");
        out2 += argv[3];
        out2 += ".txt";
        acceptfile.open(out2);



        // Monte Carlo sampling
        for (int k = 1; k <= mc_cycles; k ++){

            // Metropolis test
            Metropolis(n, E, M, w, L, accepted_configs);

            energy_sum += E;
            mag_sum += fabs(M);
            energy_sqr += E * E;
            mag_sqr += M * M;

            // write values to file
            plotfile << E << "  " << M << "  " << k << endl;
            acceptfile << accepted_configs << "  " << k << endl;

        }

        plotfile.close();
        acceptfile.close();
    }

    if (atoi(argv[4]) == 2){

        // files to write to
        ofstream probfile;

        string prob ("../files_python/probability");
        prob += argv[3];
        prob += ".txt";
        probfile.open(prob);

        // calculate energies etc. until we reach an equilibrium state
        for (int k = 1; k <= 1000; k ++){

            // Metropolis test
            Metropolis(n, E, M, w, L, accepted_configs);

        }

        // Monte Carlo sampling after equilibrium has been reached
        for (int k = 1001; k <= mc_cycles; k ++){

            // Metropolis test
            Metropolis(n, E, M, w, L, accepted_configs);

            probfile << E << " " << k << endl;

            energy_sum += E;
            energy_sqr += E * E;
        }

        probfile.close();
    }



    // write values to file for the 2x2 lattice
    if ( atoi(argv[1]) == 2 && atoi(argv[2]) >= 1000000 ){

        // mean energy
        double mean_energy = energy_sum / ((double) (mc_cycles+1));

        // mean absolute magnetic moment
        double mean_mag = mag_sum / ((double) (mc_cycles+1));

        // specific heat
        double mean_energy2 = energy_sqr / ((double) (mc_cycles+1));
        double spec_heat = ( mean_energy2 - mean_energy * mean_energy ) * beta / T_init;

        // magnetic susceptibility
        double mean_mag2 = mag_sqr / ((double) (mc_cycles+1));
        double mag_susc = ( mean_mag2 - mean_mag * mean_mag ) * beta;

        ofstream outputfile;

        outputfile.open("../files_python/2x2_output.txt");
        outputfile << "Mean energy: " << mean_energy  << endl
                   << "Mean magnetic moment: " << mean_mag<< endl
                   << "Mean energy squared: " << mean_energy2 << endl
                   << "Mean energy x mean energy: " << mean_energy* mean_energy << endl
                   << "Mean magnetic moment squared: " << mean_mag2 << endl
                   << "Mean magnetic moment x mean magnetic moment: " << mean_mag * mean_mag << endl
                   << "Specific heat: " << spec_heat << endl
                   << "Magnetic susceptibility: " << mag_susc << endl;
        outputfile.close();
    }

    return 0;
}


