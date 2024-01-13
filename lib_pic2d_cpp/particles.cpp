#include <vector>
#include <cmath>
#include <omp.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


void time_evolution_v(VectorXd& v, VectorXd& gamma, 
                      const VectorXd B_particle, const VectorXd E_particle, 
                      int n_start, int n_last, 
                      double m, double q, double dt, double c)
{
    double tmp1 = q / m * (dt / 2.0);
    double tmp4 = 1.0 / (c * c);
    

    #pragma omp parallel
    {

        Vector3d tmp_v, T, S, v_minus, v_0, v_plus;

        #pragma omp for
        for (int i = 3 * n_start; i < 3 * n_last; i+=3) {

            tmp_v = v(seq(i, i+2));

            gamma(i/3) = sqrt(1.0 + tmp_v.dot(tmp_v) * tmp4);

            T = tmp1 / gamma(i/3) * B_particle(seq(i, i+2));

            S = 2.0 / (1.0 + T.dot(T)) * T;

            v_minus = tmp_v + tmp1 * E_particle(seq(i, i+2));

            v_0 = v_minus + v_minus.cross(T);

            v_plus = v_minus + v_0.cross(S);
        
            v(seq(i, i+2)) = v_plus + tmp1 * E_particle(seq(i, i+2));
        
        }
    }
}



void time_evolution_x(VectorXd& r, 
                      const VectorXd gamma,
                      const VectorXd v,
                      int n_start, int n_last, double dt, double c)
{

    #pragma omp parallel for
    for (int i = 3 * n_start; i < 3 * n_last; i+=3){
        r(seq(i, i+2)) += v(seq(i, i+2)) / gamma(i/3) * dt;
    }
}

