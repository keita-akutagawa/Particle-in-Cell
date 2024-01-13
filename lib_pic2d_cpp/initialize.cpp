#include <vector>
#include <random>
#include <cmath>
#include <omp.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//openmpを使わないこと。STLと喧嘩しているみたい。

void set_initial_position_x(VectorXd& r, 
                          int n_start, int n_last,  
                          double x_min, double x_max, int seed)
{
    mt19937_64 mt64(seed);
    uniform_real_distribution<double> set_x(0.0, 1.0);
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        double random_value = set_x(mt64);
        if (random_value == 1.0) random_value -= 1e-10;
        r(i) = random_value * (x_max - x_min) + x_min;

    }
}

void set_initial_position_y(VectorXd& r, 
                          int n_start, int n_last,  
                          double y_min, double y_max, int seed)
{
    mt19937_64 mt64(seed);
    uniform_real_distribution<double> set_x(0.0, 1.0);
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        double random_value = set_x(mt64);
        if (random_value == 1.0) random_value -= 1e-10;
        r(i+1) = random_value * (y_max - y_min) + y_min;

    }
}

void set_initial_position_y_harris(VectorXd& r, 
                                   int n_start, int n_last,  
                                   double y_min, double y_max, double sheat_thickness, 
                                   int seed)
{
    mt19937_64 mt64(seed);
    uniform_real_distribution<double> set_x(0.0, 1.0);
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        double random_value = set_x(mt64);
        if (random_value == 1.0) random_value -= 1e-10;
        double y_position = y_max / 2.0 + sheat_thickness * atanh(2.0 * random_value - 1.0);
        if (y_position > y_max) y_position = y_max / 2.0; 
        if (y_position < y_min) y_position = y_max / 2.0; 
        r(i+1) = y_position;
    }
}

void set_initial_position_y_background(VectorXd& r, 
                                       int n_start, int n_last,  
                                       double y_min, double y_max, double sheat_thickness, 
                                       int seed)
{
    mt19937_64 mt64(seed);
    uniform_real_distribution<double> set_x(0.0, 1.0);
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        while (true)
        {
            double random_value = set_x(mt64);
            if (random_value == 1.0) random_value -= 1e-10;
            double y_position = random_value * (y_max - y_min);
            double rand_pn = set_x(mt64);
            if (rand_pn < (1.0 - 1.0 / cosh((y_position - y_max/2.0)/sheat_thickness))) {
                r(i+1) = y_position;
                break;
            }
        }
        
    }
}

void set_initial_position_z(VectorXd& r, 
                            int n_start, int n_last)
{
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        r(i+2) = 0.0;
    }
}


void set_initial_velocity(VectorXd& v,
                          int n_start, int n_last,  
                          double vx, double vy, double vz, double v_th_x, double v_th_y, double v_th_z, 
                          double c, int seed) 
{
    mt19937_64 mt64_x(seed), mt64_y(seed+1), mt64_z(seed+2);
    normal_distribution<double> set_vx(vx, v_th_x);
    normal_distribution<double> set_vy(vy, v_th_y);
    normal_distribution<double> set_vz(vz, v_th_z);
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        while (true)
        {
            double random_value = set_vx(mt64_x);
            v(i) = random_value;
            random_value = set_vy(mt64_y);
            v(i+1) = random_value;
            random_value = set_vz(mt64_z);
            v(i+2) = random_value;
            if (v(i)*v(i) + v(i+1)*v(i+1) + v(i+2)*v(i+2) >= c*c) {
                //cout << "exceeded light speed...";
                continue;
            };
            break;
        }
        
    }
}


void set_initial_velocity_x(VectorXd& v, 
                          int n_start, int n_last,  
                          double v_species, double v_th, int seed) 
{
    mt19937_64 mt64(seed);
    normal_distribution<double> set_v_comp(v_species, v_th);
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        double random_value = set_v_comp(mt64);
        v[i] = random_value;
    }
}


void set_initial_velocity_y(vector<double>& v,
                          int n_start, int n_last,  
                          double v_species, double v_th, int seed) 
{
    mt19937_64 mt64(seed);
    normal_distribution<double> set_v_comp(v_species, v_th);
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        double random_value = set_v_comp(mt64);
        v[i+1] = random_value;
    }
}


void set_initial_velocity_z(vector<double>& v,
                          int n_start, int n_last,  
                          double v_species, double v_th, int seed) 
{
    mt19937_64 mt64(seed);
    normal_distribution<double> set_v_comp(v_species, v_th);
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        double random_value = set_v_comp(mt64);
        v[i+2] = random_value;
    }
}


