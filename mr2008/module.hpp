#include <iostream>
#include <iomanip>
#include <random>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <omp.h>
//#include <filesystem>
#include <numeric>
#include <chrono>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;



void set_initial_position_x(VectorXd&, int, int, double, double, int);
void set_initial_position_y(VectorXd&, int, int, double, double, int);
void set_initial_position_y_harris(VectorXd&, int, int, double, double, double, int);
void set_initial_position_y_background(VectorXd& r, int, int, double, double, double, int);
void set_initial_position_z(VectorXd&, int, int);
void set_initial_velocity(VectorXd&, int, int, 
                          double, double, double, double, double, double,
                          double, int); 

void get_rho(VectorXd&, 
             const VectorXd, 
             int, int, double, 
             int, int, double, double);
void sort_particles(VectorXi&, 
                    VectorXd&, VectorXd&, 
                    VectorXd&, 
                    int, int, int, int);

void time_evolution_B(VectorXd&, VectorXd&, VectorXd&, 
                      const VectorXd, const VectorXd, const VectorXd,  
                      int, int, double, double, double);
void time_evolution_E(VectorXd&, VectorXd&, VectorXd&,  
                      const VectorXd, const VectorXd, const VectorXd, 
                      const VectorXd, const VectorXd, const VectorXd, 
                      int, int, double, double, double, 
                      double, double); 
void get_current_density(VectorXd&, VectorXd&, VectorXd&, 
                         VectorXd&, 
                         const VectorXd, const VectorXd,
                         int, int, 
                         double, int, int, double, double, double);
void get_particle_field(VectorXd&, VectorXd&, 
                        const VectorXd, const VectorXd, const VectorXd, 
                        const VectorXd, const VectorXd, const VectorXd, 
                        const VectorXd,  
                        int, int, int, int, double, double);
void time_evolution_v(VectorXd&, VectorXd&, 
                      const VectorXd, const VectorXd, 
                      int, int, double, double, double, double);
void time_evolution_x(VectorXd&, const VectorXd, 
                      const VectorXd, 
                      int, int, double, double);
void refrective_boudary_condition_x_left(VectorXd&, VectorXd&, int, int, double);
void refrective_boudary_condition_x_right(VectorXd&, VectorXd&, int, int, double);
void refrective_boudary_condition_y(VectorXd&, VectorXd&, int, int, double, double);
void boundary_B(VectorXd&, VectorXd&, VectorXd&, 
                int, int, double, double, double);
void boundary_E(VectorXd&, VectorXd&, VectorXd&, 
                const VectorXd, 
                int, int, double, double, 
                double, double, double);
void filter_E(VectorXd&, VectorXd&, VectorXd&, 
              const VectorXd, VectorXd&, 
              int, int, double, double, double, 
              double, double);


void io_xvKE(const VectorXd, const VectorXd, const VectorXd, 
             string, string, int, 
             int, int, 
             int, int, double, double);
void io_EBmoment(const VectorXd, const VectorXd, const VectorXd, 
                 const VectorXd, const VectorXd, const VectorXd, 
                 const VectorXd, 
                 VectorXd&, VectorXd&, 
                 VectorXd&, VectorXd&, 
                 VectorXd&, VectorXd&, 
                 string, string, int, 
                 int, int, double, double, double, double);
void get_moment(VectorXd&, VectorXd&, VectorXd&,  
                VectorXd&,
                const VectorXd, const VectorXd, 
                int, int, 
                int, int, double, double, double);