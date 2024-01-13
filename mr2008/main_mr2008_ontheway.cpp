#include <iostream>
#include <iomanip>
#include <random>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <omp.h>
//#include <filesystem>
#include <numeric>
#include <chrono>
#include "module.hpp"
#include <limits>

using namespace std;


///////////////////////////////////////////////////////////////////////


int main()
{
    string dirname = "mr2008", filename = "mr2008";
    //if(!filesystem::exists(dirname)) {
    //    filesystem::create_directory(dirname);
    //}
    double inf = numeric_limits<double>::infinity();
    int inf_int = numeric_limits<int>::infinity();

    const double c = 0.5, epsilon0 = 1.0, mu_0 = 1.0 / (epsilon0 * pow(c, 2));
    int n_e = 25, n_i = 25, step;
    int n_x, n_y;
    int n_electron, n_electron_background, n_ion, n_ion_background, n_particle;
    double dx = 1.0, dy = 1.0, x_min, x_max, y_min, y_max, dt;
    int buffer_ion, buffer_electron;
    double d; //補正用。A.B.Langdon(1992)

    double B0, B0_g;
    double t_r, T_e, T_i;
    double m_unit, r_m, m_electron, m_ion;
    double q_unit, r_q, q_ion, q_electron;
    double omega_pe, omega_pi, omega_ce, omega_ci;
    double V_A, C_S, debye_length;
    double ion_inertial_length, sheat_thickness;
    vector<double> v_ion(3, 0), v_electron(3, 0);
    double v_thi, v_the, vb_thi, vb_the;

    m_unit = 1.0;
    r_m = 1.0/25.0;
    t_r = 1.0;
    m_electron = 1 * m_unit;
    m_ion = m_electron / r_m;
    r_q = 1.0;
    B0 = sqrt(static_cast<double>(n_e)) / 1.0;
    T_i = (B0 * B0 / 2.0 / mu_0) / (n_i + n_e * t_r);
    T_e = T_i * t_r;
    debye_length = 1.0;
    q_unit = sqrt(epsilon0 * T_e / static_cast<double>(n_e)) / debye_length;
    q_electron = -1 * q_unit;
    q_ion = r_q * q_unit;
    omega_pe = sqrt(static_cast<double>(n_e) * pow(q_electron, 2) / m_electron / epsilon0);
    omega_pi = sqrt(static_cast<double>(n_i) * pow(q_ion, 2) / m_ion / epsilon0);
    omega_ce = q_electron * B0 / m_electron;
    omega_ci = q_ion * B0 / m_ion;
    V_A = B0 / sqrt(mu_0 * (static_cast<double>(n_i) * m_ion));
    C_S = sqrt(r_m * T_e);
    //debye_length = sqrt(epsilon0 * T_e / static_cast<double>(n_e) / pow(q_electron, 2));
    ion_inertial_length = c / omega_pi;
    sheat_thickness = 3.0 * ion_inertial_length;

    dt = 0.5;
    step = int(500/omega_ci/dt);
    dx = 1.0;
    dy = 1.0;
    n_x = int(ion_inertial_length * 300);
    n_y = int(ion_inertial_length * 150);
    x_min = 0.0 * dx;
    x_max = n_x * dx;
    y_min = 0.0 * dy;
    y_max = n_y * dy;
    v_thi = sqrt(T_i / m_ion);
    v_the = sqrt(T_e / m_electron);
    vb_thi = sqrt(T_i / 5.0 / m_ion);
    vb_the = sqrt(T_e / 5.0 / m_electron);
    v_electron[0] = 0.0;
    v_electron[1] = 0.0;
    v_electron[2] = c * debye_length / sheat_thickness * sqrt(2 / (1.0 + 1.0/t_r));
    v_ion[0] = -v_electron[0] / t_r;
    v_ion[1] = -v_electron[1] / t_r;
    v_ion[2] = -v_electron[2] / t_r;

    n_ion = int(n_x * n_i * 2.0 * sheat_thickness);
    n_electron = int(n_ion * r_q);
    n_ion_background = int(n_x * 0.2 * n_i * (y_max - 2.0 * sheat_thickness));
    n_electron_background = int(n_x * 0.2 * n_e * (y_max - 2.0 * sheat_thickness));
    n_particle = n_ion + n_ion_background + n_electron + n_electron_background;

    d = 1.0 / 8.0 / dt;

    cout << "Total number of particles is " << n_particle << ".\n";
    cout << "Total step is " << step <<  " (150 / Omega_ci).\n";
    
///////////////////////////////////////////////////////////////////////

    VectorXd Bx = VectorXd::Zero(n_x * n_y);
    VectorXd By = VectorXd::Zero(n_x * n_y);
    VectorXd Bz = VectorXd::Zero(n_x * n_y);
    VectorXd Bx_tmp = VectorXd::Zero(n_x * n_y);
    VectorXd By_tmp = VectorXd::Zero(n_x * n_y);
    VectorXd Bz_tmp = VectorXd::Zero(n_x * n_y);
    VectorXd Ex = VectorXd::Zero(n_x * n_y);
    VectorXd Ey = VectorXd::Zero(n_x * n_y);
    VectorXd Ez = VectorXd::Zero(n_x * n_y);
    VectorXd Ex_tmp = VectorXd::Zero(n_x * n_y);
    VectorXd Ey_tmp = VectorXd::Zero(n_x * n_y);
    VectorXd Ez_tmp = VectorXd::Zero(n_x * n_y);
    VectorXd Jx = VectorXd::Zero(n_x * n_y);
    VectorXd Jy = VectorXd::Zero(n_x * n_y);
    VectorXd Jz = VectorXd::Zero(n_x * n_y);
    VectorXd Jx_tmp = VectorXd::Zero(n_x * n_y);
    VectorXd Jy_tmp = VectorXd::Zero(n_x * n_y);
    VectorXd Jz_tmp = VectorXd::Zero(n_x * n_y);
    VectorXd rho = VectorXd::Zero(n_x * n_y);
    VectorXd phi = VectorXd::Zero(n_x * n_y);
    VectorXd F = VectorXd::Zero(n_x * n_y);
    VectorXd zeroth_moment_ion = VectorXd::Zero(n_x * n_y);
    VectorXd zeroth_moment_electron = VectorXd::Zero(n_x * n_y);
    VectorXd first_moment_ion = VectorXd::Zero(3 * n_x * n_y);
    VectorXd first_moment_electron = VectorXd::Zero(3 * n_x * n_y);
    VectorXd second_moment_ion = VectorXd::Zero(9 * n_x * n_y);
    VectorXd second_moment_electron = VectorXd::Zero(9 * n_x * n_y);
    VectorXd r = VectorXd::Zero(3 * n_particle);
    VectorXd v = VectorXd::Zero(3 * n_particle);
    VectorXd gamma = VectorXd::Zero(n_particle);
    VectorXd B_particle = VectorXd::Zero(3 * n_particle);
    VectorXd E_particle = VectorXd::Zero(3 * n_particle);
    VectorXi index_ion = VectorXi::Zero(n_particle/2);
    VectorXi index_electron = VectorXi::Zero(n_particle/2);
    VectorXd tmp_copy_double = VectorXd::Zero(n_particle/2);

    ///////////////////////////////////////////////////////////////////////

    ifstream electric_field("./mr2008/mr2008_E_40000.csv");
    ifstream magnetic_field("./mr2008/mr2008_B_40000.csv");
    ifstream position("./mr2008/mr2008_x_40000.csv");
    ifstream velocity("./mr2008/mr2008_v_40000.csv");
    ifstream ion_1_moment("./mr2008/mr2008_first_moment_ion_40000.csv");
    ifstream electron_1_moment("./mr2008/mr2008_first_moment_electron_40000.csv");
    
    string record;
    record = "step.txt";
    ofstream file_record(record);
    file_record.open(record, ios::app);
    file_record << int(0) << "step done..." << "\n";
    file_record.close();

    string str_buffer, str_conma_buffer;
    int index_i;
    
    index_i = 0;
    while (getline(position, str_buffer)) {
        istringstream i_stream(str_buffer);
        while (getline(i_stream, str_conma_buffer, ',')) {
            r[index_i] = stod(str_conma_buffer);
            index_i++;
        }
    }

    index_i = 0;
    while (getline(velocity, str_buffer)) {
        istringstream i_stream(str_buffer);
        while (getline(i_stream, str_conma_buffer, ',')) {
            v[index_i] = stod(str_conma_buffer);
            index_i++;
        }
    }

    int tmp_index = 0;
    double tmp[6];
    int index;
    while (getline(electric_field, str_buffer)) {
        istringstream i_stream(str_buffer);
        tmp_index = 0;
        while (getline(i_stream, str_conma_buffer, ',')) {
            if (tmp_index < 3) {
                tmp[tmp_index] = stod(str_conma_buffer);
            } else {
                tmp[tmp_index] = stoi(str_conma_buffer);
            };
            tmp_index++;
        }
        index = tmp[4] + tmp[3] * n_y;
        Ex(index) = tmp[0];
        Ey(index) = tmp[1];
        Ez(index) = tmp[2];
    }
    
    while (getline(magnetic_field, str_buffer)) {
        istringstream i_stream(str_buffer);
        tmp_index = 0;
        while (getline(i_stream, str_conma_buffer, ',')) {
            if (tmp_index < 3) {
                tmp[tmp_index] = stod(str_conma_buffer);
            } else {
                tmp[tmp_index] = stoi(str_conma_buffer);
            };
            tmp_index++;
        }
        index = tmp[4] + tmp[3] * n_y;
        Bx(index) = tmp[0];
        By(index) = tmp[1];
        Bz(index) = tmp[2];
    }


    double tmp4 = 1.0 / (c * c);
    for (int i = 0; i < 3 * n_particle; i+=3) {
        gamma[i/3] = sqrt(1.0 + (v[i] * v[i] + v[i+1] * v[i+1] + v[i+2] * v[i+2]) * tmp4);
    }

    while (getline(ion_1_moment, str_buffer)) {
        istringstream i_stream(str_buffer);
        tmp_index = 0;
        while (getline(i_stream, str_conma_buffer, ',')) {
            if (tmp_index < 3) {
                tmp[tmp_index] = stod(str_conma_buffer);
            } else {
                tmp[tmp_index] = stoi(str_conma_buffer);
            };
            tmp_index++;
        }
        index = tmp[4] + tmp[3] * n_y;
        first_moment_ion(index) = tmp[0];
        first_moment_ion(index + n_x*n_y) = tmp[1];
        first_moment_ion(index + 2*n_x*n_y) = tmp[2];
    }
    while (getline(electron_1_moment, str_buffer)) {
        istringstream i_stream(str_buffer);
        tmp_index = 0;
        while (getline(i_stream, str_conma_buffer, ',')) {
            if (tmp_index < 3) {
                tmp[tmp_index] = stod(str_conma_buffer);
            } else {
                tmp[tmp_index] = stoi(str_conma_buffer);
            };
            tmp_index++;
        }
        index = tmp[4] + tmp[3] * n_y;
        first_moment_electron(index) = tmp[0];
        first_moment_electron(index + n_x*n_y) = tmp[1];
        first_moment_electron(index + 2*n_x*n_y) = tmp[2];
    }

    for (int i = 0; i < n_x-1; i++) {
        for (int j = 0; j < n_y; j++) {
            Jx_tmp(j + n_y * i) = q_ion * first_moment_ion(j + n_y * i) + q_electron * first_moment_electron(j + n_y * i);
            Jy_tmp(j + n_y * i) = q_ion * first_moment_ion(j + n_y * i + n_x*n_y) + q_electron * first_moment_electron(j + n_y * i + n_x*n_y);
            Jz_tmp(j + n_y * i) = q_ion * first_moment_ion(j + n_y * i + 2*n_x*n_y) + q_electron * first_moment_electron(j + n_y * i + 2*n_x*n_y);
        }
    }

    //半整数格子点上に再定義
    for (int i = 0; i < n_x-1; i++) {
        for (int j = 0; j < n_y-1; j++) {
            Jx(j + n_y * i) = (Jx_tmp(j + n_y * i) + Jx_tmp(j + n_y * (i+1))) / 2.0;
            Jy(j + n_y * i) = (Jy_tmp(j + n_y * i) + Jy_tmp((j+1) + n_y * i)) / 2.0;
            Jz(j + n_y * i) = Jz_tmp(j + n_y * i);
        }
    }
    //とりあえずn_x-1で代用
    for (int j = 0; j < n_y; j++) {
        Jx(j + n_y * (n_x-1)) = Jx_tmp(j + n_y * (n_x-2));
        Jy(j + n_y * (n_x-1)) = Jy_tmp(j + n_y * (n_x-2));
        Jz(j + n_y * (n_x-1)) = Jz_tmp(j + n_y * (n_x-2));
    }


    for (int k = 40001; k < step+1; k++) {

        if (k % 10 == 0) {
            file_record.open(record, ios::app);
            file_record << int(k*dt) << "step done..." << "\n";
            file_record.close();
            cout << setw(5) << int(k*dt) << "step done..." << "\n"; 
        }

        //OUTPUT
        if (k % (5 * int(2/omega_ci/dt)) == 0) {
            io_xvKE(r, v, gamma, 
                    dirname, filename, k, 
                    n_ion+n_ion_background, n_electron+n_electron_background,  
                    n_x, n_y, m_ion, m_electron);
        }   

        if (k % int(2/omega_ci/dt) == 0) {
            zeroth_moment_ion = VectorXd::Zero(zeroth_moment_ion.size());
            zeroth_moment_electron = VectorXd::Zero(zeroth_moment_electron.size());
            first_moment_ion = VectorXd::Zero(first_moment_ion.size());
            first_moment_electron = VectorXd::Zero(first_moment_electron.size());
            second_moment_ion = VectorXd::Zero(second_moment_ion.size());
            second_moment_electron = VectorXd::Zero(second_moment_electron.size());
            get_moment(zeroth_moment_ion, first_moment_ion, second_moment_ion, 
                       gamma, r, v, 0, n_ion+n_ion_background, n_x, n_y, dx, dy, c);
            get_moment(zeroth_moment_electron, first_moment_electron, second_moment_electron, 
                       gamma, r, v, n_ion+n_ion_background, n_particle, n_x, n_y, dx, dy, c);
            rho = VectorXd::Zero(rho.size());
            get_rho(rho, r, 0, n_ion+n_ion_background, q_ion, n_x, n_y, dx, dy);
            get_rho(rho, r, n_ion+n_ion_background, n_particle, q_electron, n_x, n_y, dx, dy);
            io_EBmoment(Ex, Ey, Ez, Bx, By, Bz, rho, 
                        zeroth_moment_ion, zeroth_moment_electron, 
                        first_moment_ion, first_moment_electron, 
                        second_moment_ion, second_moment_electron, 
                        dirname, filename, k, 
                        n_x, n_y, dx, dy, epsilon0, mu_0);
        }     


        //STEP2
        //cout << "STEP2" << "\n";
        time_evolution_B(Bx, By, Bz, Ex, Ey, Ez, n_x, n_y, dx, dy, dt/2.0);
        boundary_B(Bx, By, Bz, n_x, n_y, dx, dy, B0_g);

        //STEP3
        //cout << "STEP3" << "\n";
        //電場・磁場を整数格子点上に再定義
        for (int i = 1; i < n_x; i++) {
            for (int j = 1; j < n_y; j++) {
                Bx_tmp(j + n_y * i) = (Bx(j + n_y * i) + Bx((j-1+n_y)%n_y + n_y * i)) / 2.0;
                By_tmp(j + n_y * i) = (By(j + n_y * i) + By(j + n_y * ((i-1+n_x)%n_x))) / 2.0;
                Bz_tmp(j + n_y * i) = (Bz(j + n_y * i) + Bz(j + n_y * ((i-1+n_x)%n_x))
                                    + Bz((j-1+n_y)%n_y + n_y * i) + Bz((j-1+n_y)%n_y + n_y * ((i-1+n_x)%n_x))) / 4.0;
                Ex_tmp(j + n_y * i) = (Ex(j + n_y * i) + Ex(j + n_y * ((i-1+n_x)%n_x))) / 2.0;
                Ey_tmp(j + n_y * i) = (Ey(j + n_y * i) + Ey((j-1+n_y)%n_y + n_y * i)) / 2.0;
                Ez_tmp(j + n_y * i) = Ez(j + n_y * i);
            }
        }
        for (int i = 0; i < n_x; i++) {
            Bx_tmp(0 + n_y * i) = Bx_tmp(1 + n_y * i);
            By_tmp(0 + n_y * i) = 0.0;
            Bz_tmp(0 + n_y * i) = Bz_tmp(1 + n_y * i);
        }
        for (int i = 0; i < n_x; i++) {
            Ex_tmp(0 + n_y * i) = 0.0;
            Ey_tmp(0 + n_y * i) = Ey_tmp(1 + n_y * i);
            Ez_tmp(0 + n_y * i) = 0.0;
        }
        for (int j = 0; j < n_y; j++) {
            Bx_tmp(j + n_y * 0) = Bx_tmp(j + n_y * 1);
            By_tmp(j + n_y * 0) = 0.0;
            Bz_tmp(j + n_y * 0) = B0_g;
        }
        for (int j = 0; j < n_y; j++) {
            Ex_tmp(j + n_y * 0) = 0.0;
            Ey_tmp(j + n_y * 0) = Ey_tmp(j + n_y * 1);
            Ez_tmp(j + n_y * 0) = Ez_tmp(j + n_y * 1);
        }
        B_particle = VectorXd::Zero(B_particle.size());
        E_particle = VectorXd::Zero(E_particle.size());
        get_particle_field(B_particle, E_particle, 
                           Bx_tmp, By_tmp, Bz_tmp, Ex_tmp, Ey_tmp, Ez_tmp, 
                           r, 0, n_ion+n_ion_background, n_x, n_y, dx, dy);
        get_particle_field(B_particle, E_particle,
                           Bx_tmp, By_tmp, Bz_tmp, Ex_tmp, Ey_tmp, Ez_tmp, 
                           r, n_ion+n_ion_background, n_particle, n_x, n_y, dx, dy);
        time_evolution_v(v, gamma,  
                         B_particle, E_particle, 
                         0, n_ion+n_ion_background, 
                         m_ion, q_ion, dt, c); 
        time_evolution_v(v, gamma,  
                         B_particle, E_particle, 
                         n_ion+n_ion_background, n_particle, 
                         m_electron, q_electron, dt, c); 

        //STEP4
        //cout << "STEP4" << "\n";
        time_evolution_x(r, gamma, v, 0, n_ion+n_ion_background, dt/2.0, c);
        time_evolution_x(r, gamma, v, n_ion+n_ion_background, n_particle, dt/2.0, c);
        refrective_boudary_condition_x_left(v, r, 0, n_particle, x_min);
        refrective_boudary_condition_x_right(v, r, 0, n_particle, x_max);
        refrective_boudary_condition_y(v, r, 0, n_particle, y_min, y_max);

        //STEP5
        //cout << "STEP5" << "\n";
        //0で初期化
        Jx_tmp = VectorXd::Zero(Jx_tmp.size());
        Jy_tmp = VectorXd::Zero(Jy_tmp.size());
        Jz_tmp = VectorXd::Zero(Jz_tmp.size());
        get_current_density(Jx_tmp, Jy_tmp, Jz_tmp, gamma, r, v, 
                            0, n_ion+n_ion_background, q_ion, n_x, n_y, dx, dy, c);
        get_current_density(Jx_tmp, Jy_tmp, Jz_tmp, gamma, r, v, 
                            n_ion+n_ion_background, n_particle, q_electron, n_x, n_y, dx, dy, c);
        //半整数格子点上に再定義
        for (int i = 0; i < n_x-1; i++) {
            for (int j = 0; j < n_y-1; j++) {
                Jx(j + n_y * i) = (Jx_tmp(j + n_y * i) + Jx_tmp(j + n_y * (i+1))) / 2.0;
                Jy(j + n_y * i) = (Jy_tmp(j + n_y * i) + Jy_tmp((j+1) + n_y * i)) / 2.0;
                Jz(j + n_y * i) = Jz_tmp(j + n_y * i);
            }
        }
        //とりあえずn_x-1で代用
        for (int j = 0; j < n_y; j++) {
            Jx(j + n_y * (n_x-1)) = Jx_tmp(j + n_y * (n_x-2));
            Jy(j + n_y * (n_x-1)) = Jy_tmp(j + n_y * (n_x-2));
            Jz(j + n_y * (n_x-1)) = Jz_tmp(j + n_y * (n_x-2));
        }


        //STEP6
        //cout << "STEP6" << "\n";
        time_evolution_B(Bx, By, Bz, Ex, Ey, Ez, n_x, n_y, dx, dy, dt/2.0);
        boundary_B(Bx, By, Bz, n_x, n_y, dx, dy, B0_g);

        
        //STEP7
        //cout << "STEP7" << "\n";
        time_evolution_x(r, gamma, v, 0, n_ion+n_ion_background, dt/2.0, c);
        time_evolution_x(r, gamma, v, n_ion+n_ion_background, n_particle, dt/2.0, c);
        refrective_boudary_condition_x_left(v, r, 0, n_particle, x_min);
        refrective_boudary_condition_x_right(v, r, 0, n_particle, x_max);
        refrective_boudary_condition_y(v, r, 0, n_particle, y_min, y_max);

        //STEP8
        //cout << "STEP8" << "\n";
        time_evolution_E(Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, n_x, n_y, dx, dy, dt, c, epsilon0);
        boundary_E(Ex, Ey, Ez, rho, n_x, n_y, dx, dy, c, epsilon0, B0_g);
        rho = VectorXd::Zero(rho.size());
        get_rho(rho, r, 0, n_ion+n_ion_background, q_ion, n_x, n_y, dx, dy);
        get_rho(rho, r, n_ion+n_ion_background, n_particle, q_electron, n_x, n_y, dx, dy);
        filter_E(Ex, Ey, Ez, rho, F, n_x, n_y, dx, dy, dt, epsilon0, d);

        if (k % 50 == 0) {
            sort_particles(index_ion, r, v, 
                           tmp_copy_double, 0, n_ion+n_ion_background, n_x, n_y);
            sort_particles(index_electron, r, v, 
                           tmp_copy_double, n_ion+n_ion_background, n_particle, n_x, n_y);
        }

    }
}
