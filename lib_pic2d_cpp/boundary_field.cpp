#include <vector>
#include <cmath>
#include <omp.h>
#include <limits>
#include <utility>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


void boundary_B(VectorXd& Bx, VectorXd& By, VectorXd& Bz, 
                int n_x, int n_y, double dx, double dy, double B0_g)
{
    for (int i = 0; i < n_x; i++) {
        By(0 + n_y * i) = 0.0;
        By(1 + n_y * i) = 0.0;
        By(n_y-1 + n_y * i) = 0.0;
        Bz(0 + n_y * i) = Bz(1 + n_y * i);
        Bz(n_y-1 + n_y * i) = Bz(n_y-2 + n_y * i);
    }
    for (int i = 0; i < n_x-1; i++) {
        Bx(0 + n_y * (i+1)) = -(By(1 + n_y * i) - By(0 + n_y * i)) / dy * dx + Bx(0 + n_y * i);
        Bx(n_y-1 + n_y * (i+1)) = -(By(n_y-1 + n_y * i) - By(n_y-2 + n_y * i)) / dy * dx + Bx(n_y-1 + n_y * i);
    }

    for (int j = 0; j < n_y; j++) {
        Bx(j + n_y * 0) = Bx(j + n_y * 1);
        Bx(j + n_y * (n_x-1)) = Bx(j + n_y * (n_x-2));
        By(j + n_y * 0) = 0.0;
        By(j + n_y * (n_x-1)) = 0.0;
        By(j + n_y * (n_x-2)) = 0.0;
        Bz(j + n_y * 0) = B0_g;
        Bz(j + n_y * (n_x-1)) = B0_g;
        Bz(j + n_y * (n_x-2)) = B0_g;
    }

}


void boundary_E(VectorXd& Ex, VectorXd& Ey, VectorXd& Ez, const VectorXd rho, 
                int n_x, int n_y, double dx, double dy, 
                double c, double epsilon0, double B0_g)
{
    for (int i = 0; i < n_x-1; i++) {
        Ex(0 + n_y * i) = 0.0;
        Ex(1 + n_y * i) = 0.0;
        Ex(n_y-1 + n_y * i) = 0.0;
        Ey(0 + n_y * i) = 0.0;
        Ey(n_y-1 + n_y * i) = 0.0;
        Ez(0 + n_y * i) = 0.0;
        Ez(1 + n_y * i) = 0.0;
        Ez(n_y-1 + n_y * i) = 0.0;
    }


    for (int j = 0; j < n_y; j++) {
        Ex(j + n_y * 0) = 0.0;
        Ex(j + n_y * (n_x-1)) = 0.0;
        Ex(j + n_y * (n_x-2)) = 0.0;
        Ey(j + n_y * 0) = Ey(j + n_y * 1);
        Ey(j + n_y * (n_x-1)) = Ey(j + n_y * (n_x-2));
        Ez(j + n_y * 0) = Ez(j + n_y * 1);
        Ez(j + n_y * (n_x-1)) = Ez(j + n_y * (n_x-2));
    }

}