#include <vector>
#include <cmath>
#include <omp.h>
#include <limits>
#include <utility>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


void filter_E(VectorXd& Ex, VectorXd& Ey, VectorXd& Ez, 
              const VectorXd rho, VectorXd& F, 
              int n_x, int n_y, double dx, double dy, double dt, 
              double epsilon0, double d)
{
    for (int i = 1; i < n_x; i++) {
        for (int j = 1; j < n_y; j++) {
            F(j + n_y * i) = ((Ex((j + n_y * i)) - Ex((j + n_y * (i-1))))/dx + (Ey((j + n_y * i)) - Ey((j-1) + n_y * i))/dy)
                    - rho(j + n_y * i) / epsilon0;
        }
    }
    for (int i = 0; i < n_x; i++) {
        F(0 + n_y * i) = 0.0;
    }
    
    
    for (int i = 0; i < n_x-1; i++) {
        for (int j = 0; j < n_y-1; j++) {
            Ex(j + n_y * i) += dt * d * (F(j + n_y * (i+1)) - F(j + n_y * i)) / dx;
            Ey(j + n_y * i) += dt * d * (F(j+1 + n_y * i) - F(j + n_y * i)) / dy;
        }
    }

    //端は0にすることでノイズを拡散させる
    for (int i = 0; i < n_x-1; i++) {
        Ex(n_y-1 + n_y * i) += dt * d * (F(n_y-1 + n_y * (i+1)) - F(n_y-1 + n_y * i)) / dx;
        Ey(n_y-1 + n_y * i) += dt * d * (0.0 - F(n_y-1 + n_y * i)) / dy;
    }
    Ex(n_y-1 + n_y * (n_x-1)) += dt * d * (0.0 - F(n_y-1 + n_y * (n_x-1)))/dx;
    Ey(n_y-1 + n_y * (n_x-1)) += dt * d * (0.0 - F(n_y-1 + n_y * (n_x-1)))/dy;

    for (int j = 0; j < n_y-1; j++) {
        Ex(j + n_y * (n_x-1)) += dt * d * (0.0 - F(j + n_y * (n_x-1))) / dx; 
        Ey(j + n_y * (n_x-1)) += dt * d * (F(j+1 + n_y * (n_x-1)) - F(j + n_y * (n_x-1))) / dy;
    }

}