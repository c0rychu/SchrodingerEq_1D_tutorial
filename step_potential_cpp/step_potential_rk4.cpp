

/****************************************************************************
  FileName     [ step_potential_rk4.cpp ]
  PackageName  [ step_potential_rk4 ]
  Synopsis     [ Solve 1-D Schrodinger equation using finite-dfference and RK4 scheme ]
  Author       [ Cory Chu ]
  Date         [ 2021 June 08 ]
  Copyright    [ Copyleft(c) 2021 GWLab, Taiwan ]
****************************************************************************/


#include <iostream>
#include <cstdio>
#include <vector>
#include <complex>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#define _USE_MATH_DEFINES
#include <cmath>


const std::complex<double> j(0.0, 1.0);
const double hbar = 1.0;
Eigen::SparseMatrix<double, Eigen::RowMajor> D2;    //RowMajor sparse matrix product can be parallelized by openmp


void D2_init(Eigen::SparseMatrix<double, Eigen::RowMajor> &d2, const unsigned int n, const double dx) {
    std::vector< Eigen::Triplet<double> > D2_tripletList;
    D2_tripletList.reserve(n*3);
    D2_tripletList.push_back(Eigen::Triplet<double>(0, 0, -2.0));
    D2_tripletList.push_back(Eigen::Triplet<double>(0, 1, 1.0));
    for (unsigned int i=1; i != n-1; ++i) {
        D2_tripletList.push_back(Eigen::Triplet<double>(i, i-1, 1.0));
        D2_tripletList.push_back(Eigen::Triplet<double>(i, i, -2.0));
        D2_tripletList.push_back(Eigen::Triplet<double>(i, i+1, 1.0));
    }
    D2_tripletList.push_back(Eigen::Triplet<double>(n-1, n-2, 1.0));
    D2_tripletList.push_back(Eigen::Triplet<double>(n-1, n-1, -2.0));

    d2.resize(n, n);
    d2.setFromTriplets(D2_tripletList.begin(), D2_tripletList.end());
    d2 /= (dx*dx);
    d2.makeCompressed();
}


// class Step_Func
/******************************************************************************/
template <class T>
class Step_Func {
    public:
        Step_Func(T x_0, T y_i, T y_f);
        Eigen::Matrix<T, Eigen::Dynamic, 1> operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x);
        T x0;
        T yi;
        T yf;
    private:
        T _step_func (const T x);
};

template <class T>
Step_Func<T>::Step_Func(T x_0, T y_i, T y_f) {
    x0 = x_0;
    yi = y_i;
    yf = y_f;
}

template <class T> 
Eigen::Matrix<T, Eigen::Dynamic, 1> Step_Func<T>::operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> y(x.size());
    for (unsigned int i=0; i != x.size(); ++i) {
        *(y.data() + i) = _step_func(*(x.data() + i));
        //y[i] = _step_func(x[i]);
    }
    return y;
}

template <class T>
T Step_Func<T>::_step_func(const T x) {
    if (x < x0) {
        return yi;
    } else {
        return yf; 
    } 
}
/******************************************************************************/


// class Gaussian_wave_packcket_Func
/******************************************************************************/
template <class T>
class Gaussian_wave_packet_Func {
    public:
        Gaussian_wave_packet_Func(T x0_, T sigma_, T kx_);
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> operator()(const Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>& x);
        T x0;
        T sigma;
        T kx;
    private:
        std::complex<T> _gaussian_wave_packet_func(const std::complex<T> x);
};

template <class T>
Gaussian_wave_packet_Func<T>::Gaussian_wave_packet_Func(T x0_, T sigma_, T kx_) {
    x0 = x0_;
    sigma = sigma_;
    kx = kx_;
}

template <class T> 
Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> Gaussian_wave_packet_Func<T>::operator()(const Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>& x) {
    Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> y(x.size());
    for (unsigned int i=0; i != x.size(); ++i) {
        *(y.data() + i) = _gaussian_wave_packet_func(*(x.data() + i));
        //y[i] = _gaussian_wave_packet_func(x[i]);
    }
    return y;
}

template <class T>
std::complex<T> Gaussian_wave_packet_Func<T>::_gaussian_wave_packet_func(const std::complex<T> x) {
    return sqrt(1.0 / (sigma * sqrt(M_PI))) * exp(-(x-x0)*(x-x0) / (2.0 * sigma * sigma)) * exp(j * kx * x);
}
/******************************************************************************/


//RHS of Schrodinger Equation
//def psi_t(t, psi):
//    return -1j * (- 0.5 * hbar / m * D2.dot(psi) + V / hbar * psi)
/******************************************************************************/
template <class T, class U>
class RHS {
    public:
        virtual T operator()(const U&, const T&) = 0;
};

class RHS_Schrodinger : public RHS< Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>, double> {
    public:
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> operator()(const double& t, const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>& psi);
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> V;
        double m;
};

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> RHS_Schrodinger::operator()(const double& t, const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>& psi) {
    return -j * (- 0.5 * hbar / m * D2 * psi + V.cwiseProduct(psi) / hbar );
}
/******************************************************************************/


// RK4 solver
/******************************************************************************/
class RK4 {
    public:
        void operator()(RHS<Eigen::VectorXcd, double>& y_t, Eigen::VectorXcd& y0);
        double ti;
        double tf;
        double dt;
        Eigen::VectorXd t;
        Eigen::MatrixXcd y;
};

void RK4::operator()(RHS<Eigen::VectorXcd, double>& y_t, Eigen::VectorXcd& y0) {
    unsigned int n_steps = floor((tf - ti) / dt);
    t.resize(n_steps + 1);
    t(0) = ti;
    y.resize(y0.size(), n_steps + 1);
    y.col(0) = y0;
    
    Eigen::VectorXcd k1(y0.size()), k2(y0.size()), k3(y0.size()), k4(y0.size()), k(y0.size());
    Eigen::VectorXcd yn(y0.size());
    double tn;

    yn = y0;
    tn = ti;
    for (unsigned int i=1; i < n_steps + 1; ++i) {
        std::cout << "RK4 iteration:" << i << std::endl;
        k1 = dt * y_t(tn, yn);
        k2 = dt * y_t(tn + dt / 2.0, yn + k1 / 2.0);
        k3 = dt * y_t(tn + dt / 2.0, yn + k2 / 2.0);
        k4 = dt * y_t(tn + dt, yn + k3);
        k = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        yn += k;
        tn += dt;
        y.col(i) = yn;
        t(i) = tn;
    }
}
/******************************************************************************/


// Export RK4 results to a file
/******************************************************************************/
class ExportRK4 {
    public:
        static void csv(FILE* f, const RK4& rk4);
        static void csv(FILE* f, const RK4& rk4, const unsigned int export_every_n_steps);
    private:
        static void _fprint_py_complex (FILE* f, const std::complex<double>& x);
};

void ExportRK4::csv(FILE* f, const RK4& rk4) {
    csv(f, rk4, 1);
}

void ExportRK4::csv(FILE* f, const RK4& rk4, const unsigned int export_every_n_steps) {
    for (int i=0; i<rk4.y.cols(); i += export_every_n_steps) {
        fprintf(f, "%f", rk4.t(i));
        fprintf(f, ", ");
        for (int k=0; k<rk4.y.rows()-1; ++k) {
            _fprint_py_complex(f, rk4.y(k, i));
            fprintf(f, ", "); 
        }
        _fprint_py_complex(f, rk4.y(rk4.y.rows()-1, i));
        fprintf(f, "\n"); 
    }
}

void ExportRK4:: _fprint_py_complex(FILE* f, const std::complex<double>& x) {
    fprintf(f, "%g%+gj", x.real(), x.imag());
}
/******************************************************************************/


// print vector
/******************************************************************************/
void fprint_vec(FILE* f, Eigen::VectorXd vec) {
    for (int i=0; i < vec.size() - 1; ++i) {
        fprintf(f, "%f, ", vec(i)); 
    }
    fprintf(f, "%f\n", vec(vec.size() - 1)); 
}
/******************************************************************************/



// main
/******************************************************************************/
int main(int argc, char **argv) {
    Eigen::initParallel();

    
    const unsigned int n_x = 1000;
    const double x_min = 0.0;
    const double x_max = 10.0;
    const double dx = (x_max - x_min) / (double)n_x;


    // D2 matrix for finite difference
    std::cout << "init D2" << std::endl;
    D2_init(D2, n_x, dx);


    // x
    Eigen::VectorXd x(n_x);
    std::cout << "init x" << std::endl;
    x.setLinSpaced(n_x, x_min, x_max - dx);


    // Step potential
    Step_Func<double> V_step(5.0, 0.0, 1000.0);
    std::cout << "init V" << std::endl;
    Eigen::VectorXd V = V_step(x);


    // Prepare gaussian wave packet
    const double x0 = 3.0; 
    const double sigma = 0.5;
    const double kx = 50.0;
    Gaussian_wave_packet_Func<double> gaussian_wave_packet(x0, sigma, kx);
    std::cout << "init psi0" << std::endl;
    Eigen::VectorXcd psi = gaussian_wave_packet(x);


    // Prepare Schrodinger equation
    RHS_Schrodinger psi_t;
    psi_t.m = 1.0;      // mass
    psi_t.V = V;


    // Prepare RK4 solver
    RK4 rk4;
    rk4.ti = 0.0;       // initial time
    rk4.tf = 0.2;       // final time
    rk4.dt = 0.0001;    // time step size


    // Run RK4 solver
    std::cout << "Start RK4" << std::endl;
    rk4(psi_t, psi);


    // Export results to a file
    FILE* oFile;
    oFile = fopen("step_potential.txt", "w");
    std::cout << "Start exporting results" << std::endl;
    fprintf(oFile, "0, ");
    fprint_vec(oFile, x);
    fprintf(oFile, "0, ");
    fprint_vec(oFile, V);
    ExportRK4::csv(oFile, rk4, 5);
    fclose(oFile);


    return 0;
}
/******************************************************************************/
