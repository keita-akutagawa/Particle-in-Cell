#include <vector>
#include <cmath>
#include <omp.h>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


void get_rho(VectorXd& rho, 
             const VectorXd r, 
             int n_start, int n_last, double q, 
             int n_x, int n_y, double dx, double dy)
{

    double cx1cy1, cx1cy2, cx2cy1, cx2cy2;
    int x_index_1, x_index_2, y_index_1, y_index_2;
    double cr_x, cr_y;

    #pragma omp parallel
    {
        
        VectorXd local_rho = VectorXd::Zero(n_x * n_y);

        #pragma omp for private(cr_x, cr_y, cx1cy1, cx1cy2, cx2cy1, cx2cy2, x_index_1, x_index_2, y_index_1, y_index_2)
        for (int i = 3 * n_start; i < 3 * n_last; i+=3) {

            x_index_1 = floor(r(i) / dx);
            x_index_2 = (x_index_1 + 1) % n_x;
            y_index_1 = floor(r(i+1) / dy);
            y_index_2 = (y_index_1 + 1) % n_y;

            cr_x = r(i) - x_index_1 * dx;
            cr_y = r(i+1) - y_index_1 * dy;

            cx1cy1 = cr_x * cr_y;
            cx1cy2 = cr_x * (1.0 - cr_y);
            cx2cy1 = (1.0 - cr_x) * cr_y;
            cx2cy2 = (1.0 - cr_x) * (1.0 - cr_y);

            local_rho(y_index_1 + n_y * x_index_1) += q * cx2cy2;
            local_rho(y_index_2 + n_y * x_index_1) += q * cx2cy1 * min({1, y_index_2});
            local_rho(y_index_1 + n_y * x_index_2) += q * cx1cy2 * min({1, x_index_2});
            local_rho(y_index_2 + n_y * x_index_2) += q * cx1cy1 * min({1, y_index_2, x_index_2});

        }

        #pragma omp critical
        {
            rho += local_rho;
        }
    }
}


/*
void get_rho_open(vector<vector<double>>& rho, 
             vector<int>& cross_r_index, 
             vector<double>& cross_cr,  
             const vector<double> cross_r,
             int n_start, int n_last, double q, 
             int n_x, int n_y, double dx, double dy)
{

    double cx1cy1, cx1cy2, cx2cy1, cx2cy2;
    int x_index_1, x_index_2, y_index_1, y_index_2;

    #pragma omp parallel for
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        cross_r_index[i] = floor(cross_r[i] / dx);
        cross_r_index[i+1] = floor(cross_r[i+1] / dy);
        cross_cr[i] = cross_r[i] / dx - cross_r_index[i];
        cross_cr[i+1] = cross_r[i+1] / dy - cross_r_index[i+1];
    }

    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {

        cx1cy1 = cross_cr[i] * cross_cr[i+1];
        cx1cy2 = cross_cr[i] * (1.0 - cross_cr[i+1]);
        cx2cy1 = (1.0 - cross_cr[i]) * cross_cr[i+1];
        cx2cy2 = (1.0 - cross_cr[i]) * (1.0 - cross_cr[i+1]);

        y_index_1 = cross_r_index[i+1];
        y_index_2 = y_index_1 + 1;

        if (y_index_2 != n_y) {

            rho[n_x-1][y_index_1] += q * cx2cy2;
            rho[n_x-1][y_index_2] += q * cx2cy1;

        } else {

            rho[n_x-1][y_index_1] += q * cx2cy2;
            rho[n_x-1][0] += q * cx2cy1;

        }
    }
}
*/

