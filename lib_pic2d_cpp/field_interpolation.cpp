#include <vector>
#include <cmath>
#include <omp.h>
#include <iostream>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


void get_particle_field(VectorXd& B_particle, VectorXd& E_particle, 
                        const VectorXd Bx_tmp, const VectorXd By_tmp, const VectorXd Bz_tmp, 
                        const VectorXd Ex_tmp, const VectorXd Ey_tmp, const VectorXd Ez_tmp, 
                        const VectorXd r, 
                        int n_start, int n_last, 
                        int n_x, int n_y, double dx, double dy)
{

    double cx1cy1, cx1cy2, cx2cy1, cx2cy2;
    int x_index_1, x_index_2, y_index_1, y_index_2;
    double cr_x, cr_y;

    #pragma omp parallel for private(cr_x, cr_y, cx1cy1, cx1cy2, cx2cy1, cx2cy2, x_index_1, x_index_2, y_index_1, y_index_2)
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

        B_particle(i) += Bx_tmp(y_index_1 + n_y * x_index_1) * cx2cy2;
        B_particle(i) += Bx_tmp(y_index_2 + n_y * x_index_1) * cx2cy1 * min({1, y_index_2});
        B_particle(i) += Bx_tmp(y_index_1 + n_y * x_index_2) * cx1cy2 * min({1, x_index_2});
        B_particle(i) += Bx_tmp(y_index_2 + n_y * x_index_2) * cx1cy1 * min({1, y_index_2, x_index_2});

        B_particle(i+1) += By_tmp(y_index_1 + n_y * x_index_1) * cx2cy2;
        B_particle(i+1) += By_tmp(y_index_2 + n_y * x_index_1) * cx2cy1 * min({1, y_index_2});
        B_particle(i+1) += By_tmp(y_index_1 + n_y * x_index_2) * cx1cy2 * min({1, x_index_2});
        B_particle(i+1) += By_tmp(y_index_2 + n_y * x_index_2) * cx1cy1 * min({1, y_index_2, x_index_2});

        B_particle(i+2) += Bz_tmp(y_index_1 + n_y * x_index_1) * cx2cy2;
        B_particle(i+2) += Bz_tmp(y_index_2 + n_y * x_index_1) * cx2cy1 * min({1, y_index_2});
        B_particle(i+2) += Bz_tmp(y_index_1 + n_y * x_index_2) * cx1cy2 * min({1, x_index_2});
        B_particle(i+2) += Bz_tmp(y_index_2 + n_y * x_index_2) * cx1cy1 * min({1, y_index_2, x_index_2});

        E_particle(i) += Ex_tmp(y_index_1 + n_y * x_index_1) * cx2cy2;
        E_particle(i) += Ex_tmp(y_index_2 + n_y * x_index_1) * cx2cy1 * min({1, y_index_2});
        E_particle(i) += Ex_tmp(y_index_1 + n_y * x_index_2) * cx1cy2 * min({1, x_index_2});
        E_particle(i) += Ex_tmp(y_index_2 + n_y * x_index_2) * cx1cy1 * min({1, y_index_2, x_index_2});

        E_particle(i+1) += Ey_tmp(y_index_1 + n_y * x_index_1) * cx2cy2;
        E_particle(i+1) += Ey_tmp(y_index_2 + n_y * x_index_1) * cx2cy1 * min({1, y_index_2});
        E_particle(i+1) += Ey_tmp(y_index_1 + n_y * x_index_2) * cx1cy2 * min({1, x_index_2});
        E_particle(i+1) += Ey_tmp(y_index_2 + n_y * x_index_2) * cx1cy1 * min({1, y_index_2, x_index_2});

        E_particle(i+2) += Ez_tmp(y_index_1 + n_y * x_index_1) * cx2cy2;
        E_particle(i+2) += Ez_tmp(y_index_2 + n_y * x_index_1) * cx2cy1 * min({1, y_index_2});
        E_particle(i+2) += Ez_tmp(y_index_1 + n_y * x_index_2) * cx1cy2 * min({1, x_index_2});
        E_particle(i+2) += Ez_tmp(y_index_2 + n_y * x_index_2) * cx1cy1 * min({1, y_index_2, x_index_2});

    }
}

