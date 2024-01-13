#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


void io_xvKE(const VectorXd r, const VectorXd v, const VectorXd gamma, 
             string dirname, string filename, int k, 
             int n_ion, int n_electron,
             int n_x, int n_y, double m_ion, double m_electron)
{
    //ios::sync_with_stdio(false);
    //cin.tie(nullptr);

    string path_x, path_v, path_ion_number, path_electron_number, path_KE, path_energy_flux;

    path_x = "./" + dirname + "/" + filename + "_x_" + to_string(k) + ".csv";
    path_v = "./" + dirname + "/" + filename + "_v_" + to_string(k) + ".csv";
    path_KE = "./" + dirname + "/" + filename + "_KE_" + to_string(k) + ".csv";
    ofstream file_x(path_x);
    ofstream file_v(path_v);
    ofstream file_KE(path_KE);
    
    double KE;
    for (int i = 0; i < 3 * n_ion; i+=3) {
        file_x << setprecision(16) << r(i) << ',' << r(i+1) << ',' << r(i+2) << "\n";
        file_v << setprecision(16) << v(i) << ',' << v(i+1) << ',' << v(i+2) << "\n";
    }
    for (int i = 3 * n_ion; i < 3 * (n_ion + n_electron); i+=3) {
        file_x << setprecision(16) << r(i) << ',' << r(i+1) << ',' << r(i+2) << "\n";
        file_v << setprecision(16) << v(i) << ',' << v(i+1) << ',' << v(i+2) << "\n";
    }

    KE = 0.0;
    KE += 1.0/2.0 * m_ion * v(seq(0, 3*n_ion-1)).dot(v(seq(0, 3*n_ion-1)));
    KE += 1.0/2.0 * m_electron * v(seq(3*n_ion, 3*(n_ion+n_electron)-1)).dot(v(seq(3*n_ion, 3*(n_ion+n_electron)-1)));
    file_KE << setprecision(16) << KE;
}


void io_EBmoment(const VectorXd Ex, const VectorXd Ey, const VectorXd Ez, 
                 const VectorXd Bx, const VectorXd By, const VectorXd Bz, 
                 const VectorXd rho, 
                 VectorXd& zeroth_moment_ion, VectorXd& zeroth_moment_electron,  
                 VectorXd& first_moment_ion, VectorXd& first_moment_electron, 
                 VectorXd& second_moment_ion, VectorXd& second_moment_electron, 
                 string dirname, string filename, int k, 
                 int n_x, int n_y, double dx, double dy, double epsilon0, double mu_0)
{
    //ios::sync_with_stdio(false);
    //cin.tie(nullptr);

    string path_E, path_B, path_energy_E, path_energy_B, path_div_E_error, path_div_B_error; 
    string path_zeroth_momemt_ion, path_zeroth_momemt_electron; 
    string path_first_momemt_ion, path_first_momemt_electron; 
    string path_second_momemt_ion, path_second_momemt_electron; 

    path_E = "./" + dirname + "/" + filename + "_E_" + to_string(k) + ".csv";
    path_B = "./" + dirname + "/" + filename + "_B_" + to_string(k) + ".csv";
    path_energy_E = "./" + dirname + "/" + filename + "_energy_E_" + to_string(k) + ".csv";
    path_energy_B = "./" + dirname + "/" + filename + "_energy_B_" + to_string(k) + ".csv";
    path_div_E_error = "./" + dirname + "/" + filename + "_div_E_error_" + to_string(k) + ".csv";
    path_div_B_error = "./" + dirname + "/" + filename + "_div_B_error_" + to_string(k) + ".csv";
    path_zeroth_momemt_ion = "./" + dirname + "/" + filename + "_zeroth_moment_ion_" + to_string(k) + ".csv";
    path_zeroth_momemt_electron = "./" + dirname + "/" + filename + "_zeroth_moment_electron_" + to_string(k) + ".csv";
    path_first_momemt_ion = "./" + dirname + "/" + filename + "_first_moment_ion_" + to_string(k) + ".csv";
    path_first_momemt_electron = "./" + dirname + "/" + filename + "_first_moment_electron_" + to_string(k) + ".csv";
    path_second_momemt_ion = "./" + dirname + "/" + filename + "_second_moment_ion_" + to_string(k) + ".csv";
    path_second_momemt_electron = "./" + dirname + "/" + filename + "_second_moment_electron_" + to_string(k) + ".csv";
    ofstream file_E(path_E);
    ofstream file_B(path_B);
    ofstream file_energy_E(path_energy_E);
    ofstream file_energy_B(path_energy_B);
    ofstream file_div_E_error(path_div_E_error);
    ofstream file_div_B_error(path_div_B_error);
    ofstream file_zeroth_moment_ion(path_zeroth_momemt_ion);
    ofstream file_zeroth_moment_electron(path_zeroth_momemt_electron);
    ofstream file_first_moment_ion(path_first_momemt_ion);
    ofstream file_first_moment_electron(path_first_momemt_electron);
    ofstream file_second_moment_ion(path_second_momemt_ion);
    ofstream file_second_moment_electron(path_second_momemt_electron);

    double energy_E, energy_B;
    energy_E = 0.0; energy_B = 0.0;
    energy_E = 1.0/2.0 * epsilon0 * (Ex.dot(Ex) + Ey.dot(Ey) + Ez.dot(Ez));
    energy_B = 1.0/2.0 / mu_0 * (Bx.dot(Bx) + By.dot(By) + Bz.dot(Bz));
    file_energy_E << setprecision(16) << energy_E;
    file_energy_B << setprecision(16) << energy_B;

    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            file_E << setprecision(16) << Ex(j + n_y * i) << ',' << Ey(j + n_y * i) << ',' << Ez(j + n_y * i) << ',' 
            << i << ',' << j << ',' << 0 << "\n";
        
            file_B << setprecision(16) << Bx(j + n_y * i) << ',' << By(j + n_y * i) << ',' << Bz(j + n_y * i) << ','
            << i << ',' << j << ',' << 0 << "\n";

            file_zeroth_moment_ion << setprecision(16) << zeroth_moment_ion(j + n_y * i) << ','
            << i << ',' << j << ',' << 0 << "\n";

            file_zeroth_moment_electron << setprecision(16) << zeroth_moment_electron(j + n_y * i) << ','
            << i << ',' << j << ',' << 0 << "\n";

            file_first_moment_ion << setprecision(16) << first_moment_ion(j + n_y * i) << ',' << first_moment_ion(j + n_y * i + n_x*n_y) << ',' << first_moment_ion(j + n_y * i + 2*n_x*n_y) << ','
            << i << ',' << j << ',' << 0 << "\n";

            file_first_moment_electron << setprecision(16) << first_moment_electron(j + n_y * i) << ',' << first_moment_electron(j + n_y * i + n_x*n_y) << ',' << first_moment_electron(j + n_y * i + 2*n_x*n_y) << ','
            << i << ',' << j << ',' << 0 << "\n";

            file_second_moment_ion << setprecision(16) 
            << second_moment_ion(j + n_y * i) << ',' << second_moment_ion(j + n_y * i + n_x*n_y) << ',' << second_moment_ion(j + n_y * i + 2*n_x*n_y) << ','
            << second_moment_ion(j + n_y * i + 3*n_x*n_y) << ',' << second_moment_ion(j + n_y * i + 4*n_x*n_y) << ',' << second_moment_ion(j + n_y * i + 5*n_x*n_y) << ','
            << second_moment_ion(j + n_y * i + 6*n_x*n_y) << ',' << second_moment_ion(j + n_y * i + 7*n_x*n_y) << ',' << second_moment_ion(j + n_y * i + 8*n_x*n_y) << ','
            << i << ',' << j << ',' << 0 << "\n";

            file_second_moment_electron << setprecision(16) 
            << second_moment_electron(j + n_y * i) << ',' << second_moment_electron(j + n_y * i + n_x*n_y) << ',' << second_moment_electron(j + n_y * i + 2*n_x*n_y) << ','
            << second_moment_electron(j + n_y * i + 3*n_x*n_y) << ',' << second_moment_electron(j + n_y * i + 4*n_x*n_y) << ',' << second_moment_electron(j + n_y * i + 5*n_x*n_y) << ','
            << second_moment_electron(j + n_y * i + 6*n_x*n_y) << ',' << second_moment_electron(j + n_y * i + 7*n_x*n_y) << ',' << second_moment_electron(j + n_y * i + 8*n_x*n_y) << ','
            << i << ',' << j << ',' << 0 << "\n";
        }
    }

    double div_B, div_E, div_E_error;
    div_B = 0.0; div_E = 0.0; div_E_error = 0.0;
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            div_B += pow((Bx(j + n_y * ((i+1)%n_x)) - Bx(j + n_y * i)) / dx + (By((j+1)%n_y + n_y * i) - By(j + n_y * i)) / dy, 2);
        }
    }
    for (int i = 1; i < n_x; i++) {
        for (int j = 1; j < n_y; j++) {
            div_E = (Ex(j + n_y * i) - Ex(j + n_y * ((i-1+n_x)%n_x))) / dx + (Ey(j + n_y * i) - Ey((j-1+n_y)%n_y + n_y * i)) / dy;
            div_E_error += pow(div_E - rho(j + n_y * i), 2);
        }
    }
    file_div_E_error << setprecision(16) << div_E_error;
    file_div_B_error << setprecision(16) << div_B;
    
    
}


void get_moment(VectorXd& zeroth_moment, VectorXd& first_moment, VectorXd& second_moment,  
                VectorXd& gamma,
                const VectorXd r, const VectorXd v, 
                int n_start, int n_last, 
                int n_x, int n_y, double dx, double dy, double c)
{

    double tmp, tmp0, tmp1, tmp2;
    double cx1cy1, cx1cy2, cx2cy1, cx2cy2;
    int x_index_1, x_index_2, y_index_1, y_index_2;
    double cr_x, cr_y;
    
    #pragma omp parallel
    {

        VectorXd local_zeroth_moment = VectorXd::Zero(n_x * n_y);
        VectorXd local_first_moment = VectorXd::Zero(n_x * n_y * 3);
        VectorXd local_second_moment = VectorXd::Zero(n_x * n_y * 9);

        #pragma omp for private(tmp, tmp0, tmp1, tmp2, cr_x, cr_y, cx1cy1, cx1cy2, cx2cy1, cx2cy2, x_index_1, x_index_2, y_index_1, y_index_2)
        for (int i = 3 * n_start; i < 3 * n_last; i+=3) {

            tmp = 1.0 / gamma(i/3);
            tmp0 = tmp * v(i);
            tmp1 = tmp * v(i+1);
            tmp2 = tmp * v(i+2);

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


            local_zeroth_moment(y_index_1 + n_y * x_index_1) += cx2cy2;
            local_zeroth_moment(y_index_2 + n_y * x_index_1) += cx2cy1 * min({1, y_index_2});
            local_zeroth_moment(y_index_1 + n_y * x_index_2) += cx1cy2 * min({1, x_index_2});
            local_zeroth_moment(y_index_2 + n_y * x_index_2) += cx1cy1 * min({1, y_index_2, x_index_2});
            //////////
            local_first_moment(y_index_1 + n_y * x_index_1) += tmp0 * cx2cy2;
            local_first_moment(y_index_2 + n_y * x_index_1) += tmp0 * cx2cy1 * min({1, y_index_2});
            local_first_moment(y_index_1 + n_y * x_index_2) += tmp0 * cx1cy2 * min({1, x_index_2});
            local_first_moment(y_index_2 + n_y * x_index_2) += tmp0 * cx1cy1 * min({1, y_index_2, x_index_2});

            local_first_moment(y_index_1 + n_y * x_index_1 + n_x*n_y) += tmp1 * cx2cy2;
            local_first_moment(y_index_2 + n_y * x_index_1 + n_x*n_y) += tmp1 * cx2cy1 * min({1, y_index_2});
            local_first_moment(y_index_1 + n_y * x_index_2 + n_x*n_y) += tmp1 * cx1cy2 * min({1, x_index_2});
            local_first_moment(y_index_2 + n_y * x_index_2 + n_x*n_y) += tmp1 * cx1cy1 * min({1, y_index_2, x_index_2});

            local_first_moment(y_index_1 + n_y * x_index_1 + 2*n_x*n_y) += tmp2 * cx2cy2;
            local_first_moment(y_index_2 + n_y * x_index_1 + 2*n_x*n_y) += tmp2 * cx2cy1 * min({1, y_index_2});
            local_first_moment(y_index_1 + n_y * x_index_2 + 2*n_x*n_y) += tmp2 * cx1cy2 * min({1, x_index_2});
            local_first_moment(y_index_2 + n_y * x_index_2 + 2*n_x*n_y) += tmp2 * cx1cy1 * min({1, y_index_2, x_index_2});
            //////////
            local_second_moment(y_index_1 + n_y * x_index_1) += tmp0 * tmp0 * cx2cy2;
            local_second_moment(y_index_2 + n_y * x_index_1) += tmp0 * tmp0 * cx2cy1 * min({1, y_index_2});
            local_second_moment(y_index_1 + n_y * x_index_2) += tmp0 * tmp0 * cx1cy2 * min({1, x_index_2});
            local_second_moment(y_index_2 + n_y * x_index_2) += tmp0 * tmp0 * cx1cy1 * min({1, y_index_2, x_index_2});

            local_second_moment(y_index_1 + n_y * x_index_1 + n_x*n_y) += tmp0 * tmp1 * cx2cy2;
            local_second_moment(y_index_2 + n_y * x_index_1 + n_x*n_y) += tmp0 * tmp1 * cx2cy1 * min({1, y_index_2});
            local_second_moment(y_index_1 + n_y * x_index_2 + n_x*n_y) += tmp0 * tmp1 * cx1cy2 * min({1, x_index_2});
            local_second_moment(y_index_2 + n_y * x_index_2 + n_x*n_y) += tmp0 * tmp1 * cx1cy1 * min({1, y_index_2, x_index_2});

            local_second_moment(y_index_1 + n_y * x_index_1 + 2*n_x*n_y) += tmp0 * tmp2 * cx2cy2;
            local_second_moment(y_index_2 + n_y * x_index_1 + 2*n_x*n_y) += tmp0 * tmp2 * cx2cy1 * min({1, y_index_2});
            local_second_moment(y_index_1 + n_y * x_index_2 + 2*n_x*n_y) += tmp0 * tmp2 * cx1cy2 * min({1, x_index_2});
            local_second_moment(y_index_2 + n_y * x_index_2 + 2*n_x*n_y) += tmp0 * tmp2 * cx1cy1 * min({1, y_index_2, x_index_2});

            local_second_moment(y_index_1 + n_y * x_index_1 + 3*n_x*n_y) += tmp1 * tmp0 * cx2cy2;
            local_second_moment(y_index_2 + n_y * x_index_1 + 3*n_x*n_y) += tmp1 * tmp0 * cx2cy1 * min({1, y_index_2});
            local_second_moment(y_index_1 + n_y * x_index_2 + 3*n_x*n_y) += tmp1 * tmp0 * cx1cy2 * min({1, x_index_2});
            local_second_moment(y_index_2 + n_y * x_index_2 + 3*n_x*n_y) += tmp1 * tmp0 * cx1cy1 * min({1, y_index_2, x_index_2});

            local_second_moment(y_index_1 + n_y * x_index_1 + 4*n_x*n_y) += tmp1 * tmp1 * cx2cy2;
            local_second_moment(y_index_2 + n_y * x_index_1 + 4*n_x*n_y) += tmp1 * tmp1 * cx2cy1 * min({1, y_index_2});
            local_second_moment(y_index_1 + n_y * x_index_2 + 4*n_x*n_y) += tmp1 * tmp1 * cx1cy2 * min({1, x_index_2});
            local_second_moment(y_index_2 + n_y * x_index_2 + 4*n_x*n_y) += tmp1 * tmp1 * cx1cy1 * min({1, y_index_2, x_index_2});

            local_second_moment(y_index_1 + n_y * x_index_1 + 5*n_x*n_y) += tmp1 * tmp2 * cx2cy2;
            local_second_moment(y_index_2 + n_y * x_index_1 + 5*n_x*n_y) += tmp1 * tmp2 * cx2cy1 * min({1, y_index_2});
            local_second_moment(y_index_1 + n_y * x_index_2 + 5*n_x*n_y) += tmp1 * tmp2 * cx1cy2 * min({1, x_index_2});
            local_second_moment(y_index_2 + n_y * x_index_2 + 5*n_x*n_y) += tmp1 * tmp2 * cx1cy1 * min({1, y_index_2, x_index_2});

            local_second_moment(y_index_1 + n_y * x_index_1 + 6*n_x*n_y) += tmp2 * tmp0 * cx2cy2;
            local_second_moment(y_index_2 + n_y * x_index_1 + 6*n_x*n_y) += tmp2 * tmp0 * cx2cy1 * min({1, y_index_2});
            local_second_moment(y_index_1 + n_y * x_index_2 + 6*n_x*n_y) += tmp2 * tmp0 * cx1cy2 * min({1, x_index_2});
            local_second_moment(y_index_2 + n_y * x_index_2 + 6*n_x*n_y) += tmp2 * tmp0 * cx1cy1 * min({1, y_index_2, x_index_2});

            local_second_moment(y_index_1 + n_y * x_index_1 + 7*n_x*n_y) += tmp2 * tmp1 * cx2cy2;
            local_second_moment(y_index_2 + n_y * x_index_1 + 7*n_x*n_y) += tmp2 * tmp1 * cx2cy1 * min({1, y_index_2});
            local_second_moment(y_index_1 + n_y * x_index_2 + 7*n_x*n_y) += tmp2 * tmp1 * cx1cy2 * min({1, x_index_2});
            local_second_moment(y_index_2 + n_y * x_index_2 + 7*n_x*n_y) += tmp2 * tmp1 * cx1cy1 * min({1, y_index_2, x_index_2});

            local_second_moment(y_index_1 + n_y * x_index_1 + 8*n_x*n_y) += tmp2 * tmp2 * cx2cy2;
            local_second_moment(y_index_2 + n_y * x_index_1 + 8*n_x*n_y) += tmp2 * tmp2 * cx2cy1 * min({1, y_index_2});
            local_second_moment(y_index_1 + n_y * x_index_2 + 8*n_x*n_y) += tmp2 * tmp2 * cx1cy2 * min({1, x_index_2});
            local_second_moment(y_index_2 + n_y * x_index_2 + 8*n_x*n_y) += tmp2 * tmp2 * cx1cy1 * min({1, y_index_2, x_index_2});
        }

        #pragma omp critical
        {
            zeroth_moment += local_zeroth_moment;
            first_moment += local_first_moment;
            second_moment += local_second_moment;
        }
    }
}



