#ifndef FUNCTIONS_MPI_H
#define FUNCTIONS_MPI_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void CreateRandomMatrix(int n, mat &L) ;
double CalculateEnergy(int n, mat &L) ;
double CalculateMagneticMoment(mat &L) ;
void Metropolis(int n, double &E, double &M, double *w, mat &L) ;

#endif // FUNCTIONS_MPI_H
