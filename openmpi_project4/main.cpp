#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <fstream>
#include "functions_mpi.h"
#include <mpi.h>
#include <string>

using namespace std;
using namespace arma;


int main(){

    // parameters
    int n = 40;
    int mc_cycles = 100000;
    double T_init = 1.6;
    double T_final = 3.0;
    double deltaT = 0.05;

    // create different filenames based on the lattice size 'n'
    ofstream mpifile;
    string filename ("../files_python/mpi_lattice");
    filename += to_string(n);
    filename += ".txt";
    mpifile.open(filename);

    // number of processors
    int numprocs;
    // rank of calling process
    int my_rank;

    double w[17];
    double average[4];
    double average_tot[4];

    mat L = zeros<mat>(n, n);

    // make each random number be random every time we call rand()
    srand (time(NULL));

    // MPI initializations
    int initialized, finalized;

    MPI_Initialized(&initialized);
    if (!initialized){
        // initialize MPI by NULL, and not by input arguments(we don't have any)
        MPI_Init(NULL, NULL);
    }

    // returns size of number of processors
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    // returns the rank of the calling process
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


    // define where each processor begins and ends its loop
    int intervals = mc_cycles / numprocs;
    int begin_loop = my_rank * intervals + 1;
    int end_loop = ( my_rank + 1 ) * intervals;

    if ( (my_rank == numprocs - 1) && (end_loop < mc_cycles) ){
        end_loop = mc_cycles;
    }

    // broadcast variables to all nodes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T_init, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T_final, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&deltaT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    // timing variables
    double T_start;
    double T_end;
    double T_total;

    T_start = MPI_Wtime();

    for (double temp = T_init; temp <= T_final; temp+=deltaT){

        // create a random matrix
        CreateRandomMatrix(n, L);

        // initialize energy and magnetization
        double E = CalculateEnergy(n, L);
        double M = CalculateMagneticMoment(L);


        // initialize array
        for (int dE = -8; dE <= 8; dE++){ w[dE + 8] = 0; }
        // every eight element is an exponential of energy
        for (int dE = -8; dE <= 8; dE+=4){ w[dE + 8] = exp(- dE / temp); }

        // initialize arrays
        for (int i = 0; i < 4; i++){ average[i] = 0; }
        for (int i = 0; i < 4; i++){ average_tot[i] = 0; }


        for (int cycles = begin_loop; cycles <= end_loop; cycles++){

            Metropolis(n, E, M, w, L);

            average[0] += E;
            average[1] += fabs(M);
            average[2] += E * E;
            average[3] += M * M;

        }

        for (int i = 0; i < 4; i++){

            // MPI_Reduce takes an array of input element on each process
            // and returns an array of output elements to the root process
            MPI_Reduce(&average[i], &average_tot[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }

        // calculate variables
        double norm = 1.0 / ((double) mc_cycles);
        double mean_energy = average_tot[0] * norm;
        double mean_mag = average_tot[1] * norm;
        double energy_sqr = average_tot[2] * norm;
        double mag_sqr = average_tot[3] * norm;
        double spec_heat = ( energy_sqr - (mean_energy*mean_energy) ) / (temp * temp);
        double mag_susc = ( mag_sqr - (mean_mag*mean_mag) ) / (temp);

        // write to file at the end of every process
        if (my_rank == 0){
            mpifile << n << " " << temp << " " << mean_energy << " " <<
                       mean_mag<< " " << spec_heat << " " << mag_susc << endl;
        }
    }

    T_end = MPI_Wtime();
    T_total = T_end - T_start;

    if (my_rank == 0){
        cout << "Time = " << T_total << " on number of processors: " << numprocs << endl;
    }

    // finalize MPI
    MPI_Finalized(&finalized);
    if (!finalized){
        MPI_Finalize();
    }

    return 0;
}


