#include <vector>
#include <cmath>
#include <omp.h>
#include <iostream>

using namespace std;


void poisson_solver_jacobi(vector<vector<double>>& phi, 
                           const vector<vector<double>> rho, int iteration,
                           int n_x, int n_y, double dx, double dy, double epsilon0)
{
    double tmp = 1.0 / (2.0 * (1.0/pow(dx, 2) + 1.0/pow(dy, 2)));
    double tmp_x = 1.0 / pow(dx, 2);
    double tmp_y = 1.0 / pow(dy, 2);
    
    //Gauss-Seidel法
    for (int iter = 0; iter < iteration; iter++) {
        for (int i = 0; i < n_x; i++) {
            for (int j = 0; j < n_y; j++) {
                phi[i][j] = tmp * (rho[i][j]/epsilon0 
                          + (phi[(i-1+n_x)%n_x][j] + phi[(i+1)%n_x][j]) * tmp_x
                          + (phi[i][(j-1+n_y)%n_y] + phi[i][(j+1)%n_y]) * tmp_y);
            }
        }
        for (int j = 0; j < n_y; j++) {
            phi[0][j] = 0.0;
            phi[n_x-1][j] = 0.0;
        }
        for (int i = 0; i < n_x; i++) {
            phi[i][0] = phi[i][1];
            phi[i][n_y-1] = phi[i][n_y-2];
        }
    }
}


void get_E(vector<vector<vector<double>>>& E, vector<vector<double>> phi, 
           int n_x, int n_y, double dx, double dy)
{
    //端は境界条件の方で決まるので気にしなくてよい
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            E[0][i][j] = -(phi[(i+1)%n_x][j] - phi[i][j]) / dx;
            E[1][i][j] = -(phi[i][(j+1)%n_y] - phi[i][j]) / dy;
        }
    }
}



