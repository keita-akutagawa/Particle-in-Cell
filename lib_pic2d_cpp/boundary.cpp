#include <vector>
#include <cmath>
#include <omp.h>
#include <limits>
#include <utility>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;



void refrective_boudary_condition_x_left(VectorXd& v, VectorXd& r,  
                                         int n_start, int n_last, double x_min)
{
    #pragma omp parallel for
    for (int i = 3 * n_start; i < 3 * n_last; i+=3){
        if (r(i) <= x_min) {
            v(i) = -v(i);
            r(i) = 2.0 * x_min - r(i);
        }
    }
}


void refrective_boudary_condition_x_right(VectorXd& v, VectorXd& r,  
                                         int n_start, int n_last, double x_max)
{
    #pragma omp parallel for
    for (int i = 3 * n_start; i < 3 * n_last; i+=3){
        if (r(i) >= x_max) {
            v(i) = -v(i);
            r(i) = 2.0 * x_max - r(i);
        }
    }
}


void refrective_boudary_condition_y(VectorXd& v, VectorXd& r,  
                                    int n_start, int n_last, double y_min, double y_max)
{
    #pragma omp parallel for
    for (int i = 3 * n_start; i < 3 * n_last; i+=3){
        if (r(i+1) >= y_max) {
            v(i+1) = -v(i+1);
            r(i+1) = 2.0 * y_max - r(i+1);
        } else if (r(i+1) <= y_min) {
            v(i+1) = -v(i+1);
            r(i+1) = 2.0 * y_min - r(i+1);
        }
    }
}


