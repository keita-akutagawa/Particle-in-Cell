#include <vector>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


void time_evolution_B(VectorXd& Bx, VectorXd& By, VectorXd& Bz, 
                      const VectorXd Ex, const VectorXd Ey, const VectorXd Ez,  
                      int n_x, int n_y, double dx, double dy, double dt) 
{
    int i, j = 0;
    #pragma omp parallel for private(j)
    for (i = 0; i < n_x-1; i++) {
        for (j = 0; j < n_y-1; j++) {
            Bx(j + n_y * i) += -(Ez((j+1) + n_y * i) - Ez(j + n_y * i)) / dy * dt;
            By(j + n_y * i) += (Ez(j + n_y * (i+1)) - Ez(j + n_y * i)) / dx * dt;
            Bz(j + n_y * i) += (-(Ey(j + n_y * (i+1)) - Ey(j + n_y * i)) / dx
                        + (Ex((j+1) + n_y * i) - Ex(j + n_y * i)) / dy) * dt;
        }
    }
}


void time_evolution_E(VectorXd& Ex, VectorXd& Ey, VectorXd& Ez, 
                      const VectorXd Bx, const VectorXd By, const VectorXd Bz, 
                      const VectorXd Jx, const VectorXd Jy, const VectorXd Jz,  
                      int n_x, int n_y, double dx, double dy, double dt, 
                      double c, double epsilon0) 
{
    int i, j = 0;
    #pragma omp parallel for private(j)
    for (i = 1; i < n_x; i++) {
        for (j = 1; j < n_y; j++) {
            Ex(j + n_y * i) += (-Jx(j + n_y * i) / epsilon0
                        + c * c * (Bz(j + n_y * i) - Bz((j-1) + n_y * i)) / dy) * dt;
            Ey(j + n_y * i) += (-Jy(j + n_y * i) / epsilon0 
                        - c * c * (Bz(j + n_y * i) - Bz(j + n_y * (i-1))) / dx) * dt;
            Ez(j + n_y * i) += (-Jz(j + n_y * i) / epsilon0 
                        + c * c * ((By(j + n_y * i) - By(j + n_y * (i-1))) / dx
                        - (Bx(j + n_y * i) - Bx((j-1) + n_y * i)) / dy)) * dt;
        }
    }
}
