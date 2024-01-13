#include <vector>
#include <cmath>
#include <omp.h>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


void sort_particles(VectorXi& index_species,
                    VectorXd& r, VectorXd& v,   
                    VectorXd& tmp_copy_double, 
                    int n_start, int n_last, 
                    int n_x, int n_y)
{
    #pragma omp parallel for
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        tmp_copy_double(i/3 - n_start) = r(i+1) + r(i) * n_y;
    }

    #pragma omp parallel for
    for (int i = 0; i < n_last-n_start; i++) {
        index_species(i) = i;
    }
    sort(index_species.begin(), index_species.begin() + n_last - n_start, 
         [&](int a, int b){return tmp_copy_double(a) < tmp_copy_double(b);});


    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        tmp_copy_double(i/3 - n_start) = r(i);
    }

    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        r(i) = tmp_copy_double(index_species(i/3-n_start));
    }

    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {   
        tmp_copy_double(i/3 - n_start) = r(i+1);
    }

    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        r(i+1) = tmp_copy_double(index_species(i/3-n_start));
    }

    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {   
        tmp_copy_double(i/3 - n_start) = r(i+2);
    }

    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        r(i+2) = tmp_copy_double(index_species(i/3-n_start));
    }

    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {        
        tmp_copy_double(i/3 - n_start) = v(i);
    }

    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        v(i) = tmp_copy_double(index_species(i/3-n_start));
    }

    #pragma omp parallel for    
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        tmp_copy_double(i/3 - n_start) = v(i+1);
    }

    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        v(i+1) = tmp_copy_double(index_species(i/3-n_start));
    }

    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {   
        tmp_copy_double(i/3 - n_start) = v(i+2);
    }

    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        v(i+2) = tmp_copy_double(index_species(i/3-n_start));
    }

}


