#ifndef FUNCTION_H
#define FUNCTION_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void CreateRandomMatrix(int n, mat &L) ;
void CreateMatrix(int n, mat &L) ;
double CalculateEnergy(int n, mat &L) ;
double CalculateMagneticMoment(mat &L) ;
void Metropolis(int n, double &E, double &M, double *w, mat &L) ;

#endif // FUNCTION_H
