//
//  horizon-fluxes.cpp
//  This C++ program computes P_res and other quantities defined in the paper "Smooth binary evolution from wide resonances in boson clouds"
//
//  Compile this way (gsl library directory should be changed to its location on the user's machine):
//
//  g++ -Wall -I/opt/homebrew/Cellar/gsl/2.7.1/include -c XXX.cpp
//  g++ XXX.o -L/opt/homebrew/Cellar/gsl/2.7.1/lib -lgsl -lgslcblas
//  ./a.out
//
//  Created by Giovanni Maria Tomaselli on 29/05/25.
//

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>
#include <iostream>
#include <float.h>
#include <fstream>
#include <assert.h>
#include <chrono>
#include <complex>
#include <map>
#include <set>
#include <boost/multiprecision/cpp_complex.hpp>
#include "gsl/gsl_sf_coulomb.h"
#include "gsl/gsl_sf_coupling.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_roots.h"

// Units:
// M = 1 (note that this is different from https://github.com/thomasspieksma/GrAB, where we used rs = 1 instead)

using namespace std::literals::complex_literals;
using namespace boost::multiprecision;

// Change this line to use cpp_complex_100, in order increase numerical precision in the Leaver method (needed for high-\ell states, but makes the code much slower):
using Complex_high_prec = std::complex<double>;//cpp_complex_100;//

inline long double factorial (int n)
{
    if (n == 0) return 1;
    else return n * factorial(n-1);
}

double Wigner_small_d (int j, int mprime, int m, double beta)
{
    long double x = 0;
    int smin = std::max(0, m-mprime);
    int smax = std::min(j+m, j-mprime);
    
    for (int s = smin; s <= smax; s++)
        x += (pow(-1,mprime-m+s) * pow(cos(beta/2),2*j+m-mprime-2*s) * pow(sin(beta/2),mprime-m+2*s) ) / (factorial(j+m-s)*factorial(s)*factorial(mprime-m+s)*factorial(j-mprime-s));
        
    return x * sqrt(factorial(j+mprime)*factorial(j-mprime)*factorial(j+m)*factorial(j-m));
}

// Angular overlap integral
double I_Omega (int lprime, int l_star, int l, int mprime, int m_star, int m)
{
    return pow(-1,mprime+m_star) * sqrt(((2*lprime+1)*(2*l_star+1)*(2*l+1))/(4*M_PI)) * gsl_sf_coupling_3j(2*lprime, 2*l_star, 2*l, 0, 0, 0) * gsl_sf_coupling_3j(2*lprime, 2*l_star, 2*l, -2*mprime, -2*m_star, 2*m);
}

// Wavefunction of bound states, normalized such that: integral_0^\infty Rbound^2 r^2 dr = 1
double Rbound (int n, int l, double r, double alpha)
{
    double rbohr = 1 / pow(alpha,2);
            
    return gsl_sf_hydrogenicR(n,l,1/rbohr,r);
}

// Analytical approximation of the real part of the energy eigenvalue
double Energy_bound (int n, int l, int m, double alpha, double atilde)
{
    double fnl = 0; //(2-4*n)/n; // see eqs (A.18) of 1804.03208
    if (l >= 0) fnl = (double) 2/n - (double) 6/(2*l+1);
    
    double hl = 0;
    if (l > 0) hl = (double) 16/(2*l*(2*l+1)*(2*l+2));
                            
    return alpha * ( - pow(alpha/n,2)/2 - pow(alpha/n,4)/8 + fnl * pow(alpha,4) / pow(n,3) + hl*atilde*m * pow(alpha,5) / pow(n,3) );
}

// Analytical approximation of the imaginary part of the energy eigenvalue
double Gamma_Detweiler (int n, int l, int m, double alpha, double atilde)
{
    double factor1 = pow(2,4*l+1) * factorial(n+l) / (pow(n,2*l+4) * factorial(n-l-1));
    double factor2 = factorial(l) / (factorial(2*l) * factorial(2*l+1));
    double Cnl = factor1 * pow(factor2,2);
    
    double r_plus = 1 + sqrt(1-pow(atilde,2));
    double Omega_plus = atilde / (2 * r_plus);
    
    double omeganlm = alpha + Energy_bound(n,l,m,alpha,atilde);
    
    double glm = 1.0;
    for (int k = 1; k <= l; k++)
        glm *= k*k * (1-pow(atilde,2)) + pow(atilde * m - 2 * r_plus * omeganlm, 2);
        
    return 2*r_plus * Cnl * glm * (m*Omega_plus - omeganlm) * pow(alpha, 4*l+5);
}

// Taylor series approximation of the eigenvalue of the spheroidal harmonics
Complex_high_prec SpheroidalEigenvalue (int l, int m, Complex_high_prec c_squared)
{
    Complex_high_prec a0 = l*(1 + l);
    Complex_high_prec a1 = - (2*(-1 + l + pow(l,2) + pow(m,2)))/((-1 + 2*l)*(3 + 2*l));
    Complex_high_prec a2 = ((-1 + l - m)*(l - m)*(-1 + l + m)*(l + m))/(2.*(-3 + 2*l)*pow(-1 + 2*l,3)*(1 + 2*l)) - ((1 + l - m)*(2 + l - m)*(1 + l + m)*(2 + l + m))/(2.*(1 + 2*l)*pow(3 + 2*l,3)*(5 + 2*l));
    Complex_high_prec a3 = - (4*(-1 + 4*pow(m,2))*(-15 + 121*l + 334*pow(l,2) + 130*pow(l,3) - 595*pow(l,4) - 568*pow(l,5) + 184*pow(l,6) + 320*pow(l,7) + 80*pow(l,8) - 690*pow(m,2) + 274*l*pow(m,2) - 62*pow(l,2)*pow(m,2) - 896*pow(l,3)*pow(m,2) - 1008*pow(l,4)*pow(m,2) - 672*pow(l,5)*pow(m,2) - 224*pow(l,6)*pow(m,2) + 705*pow(m,4) + 1000*l*pow(m,4) + 1144*pow(l,2)*pow(m,4) + 288*pow(l,3)*pow(m,4) + 144*pow(l,4)*pow(m,4)))/((-5 + 2*l)*(-3 + 2*l)*pow(-1 + 2*l,5)*pow(3 + 2*l,5)*(5 + 2*l)*(7 + 2*l));
    Complex_high_prec a4 = (2*(110565 - 1930851*l + 28394820*pow(l,2) + 25369092*pow(l,3) - 68210954*pow(l,4) - 69765510*pow(l,5) + 38630168*pow(l,6) + 73666140*pow(l,7) + 19843829*pow(l,8) - 30142668*pow(l,9) - 25803260*pow(l,10) + 1887808*pow(l,11) + 8380320*pow(l,12) + 1208704*pow(l,13) - 1399168*pow(l,14) - 419840*pow(l,15) + 78080*pow(l,16) + 46080*pow(l,17) + 5120*pow(l,18) - 14118300*pow(m,2) - 82386180*l*pow(m,2) - 108866700*pow(l,2)*pow(m,2) + 38392596*pow(l,3)*pow(m,2) + 182354996*pow(l,4)*pow(m,2) + 50988588*pow(l,5)*pow(m,2) - 127055740*pow(l,6)*pow(m,2) + 23756352*pow(l,7)*pow(m,2) + 94855248*pow(l,8)*pow(m,2) - 61697920*pow(l,9)*pow(m,2) - 86969216*pow(l,10)*pow(m,2) + 13805568*pow(l,11)*pow(m,2) + 41018880*pow(l,12)*pow(m,2) + 10644480*pow(l,13)*pow(m,2) - 3640320*pow(l,14)*pow(m,2) - 2064384*pow(l,15)*pow(m,2) - 258048*pow(l,16)*pow(m,2) + 302556870*pow(m,4) - 80850798*l*pow(m,4) + 433926720*pow(l,2)*pow(m,4) + 586980756*pow(l,3)*pow(m,4) - 634469610*pow(l,4)*pow(m,4) - 617676120*pow(l,5)*pow(m,4) + 595410360*pow(l,6)*pow(m,4) + 595484928*pow(l,7)*pow(m,4) - 52179264*pow(l,8)*pow(m,4) - 255836928*pow(l,9)*pow(m,4) - 148985088*pow(l,10)*pow(m,4) - 22828032*pow(l,11)*pow(m,4) + 21634560*pow(l,12)*pow(m,4) + 11741184*pow(l,13)*pow(m,4) + 1677312*pow(l,14)*pow(m,4) - 783650700*pow(m,6) - 251355420*l*pow(m,6) - 181390860*pow(l,2)*pow(m,6) - 289807200*pow(l,3)*pow(m,6) - 1053358960*pow(l,4)*pow(m,6) - 623736960*pow(l,5)*pow(m,6) + 572298880*pow(l,6)*pow(m,6) + 665272320*pow(l,7)*pow(m,6) + 141258240*pow(l,8)*pow(m,6) - 48921600*pow(l,9)*pow(m,6) - 41999360*pow(l,10)*pow(m,6) - 17571840*pow(l,11)*pow(m,6) - 2928640*pow(l,12)*pow(m,6) + 495101565*pow(m,8) + 853977924*l*pow(m,8) + 815701860*pow(l,2)*pow(m,8) - 209971008*pow(l,3)*pow(m,8) - 427288928*pow(l,4)*pow(m,8) - 353777280*pow(l,5)*pow(m,8) - 58434944*pow(l,6)*pow(m,8) + 60017664*pow(l,7)*pow(m,8) + 26286336*pow(l,8)*pow(m,8) + 7521280*pow(l,9)*pow(m,8) + 1504256*pow(l,10)*pow(m,8)))/((-7 + 2*l)*(-5 + 2*l)*pow(-3 + 2*l,3)*pow(-1 + 2*l,7)*pow(3 + 2*l,7)*pow(5 + 2*l,3)*(7 + 2*l)*(9 + 2*l));
    
    return a0 + a1 * c_squared + a2 * pow(c_squared,2) + a3 * pow(c_squared,3) + a4 * pow(c_squared,4);
}

// Newton-Raphson method for real functions
double real_newton(std::function<double(double)> f, double x0, double tol = 1e-10, int min_iter = 3, int max_iter = 100, double h_eps = 1e-10)
{
    double x = x0;

    for (int i = 0; i < max_iter; ++i) {
        double fx = f(x);

        if (std::abs(fx) < tol && i > min_iter)
            return x;

        // Approximate derivative f'(x) using symmetric finite difference
        double df = (f(x + h_eps) - f(x - h_eps)) / (2.0 * h_eps);

        if (std::abs(df) < tol)
            throw std::runtime_error("Derivative too small; possible division by zero.");

        x = x - fx / df;
    }

    throw std::runtime_error("Max iterations reached without convergence.");
}

// Newton-Raphson method for complex functions
Complex_high_prec complex_newton(std::function<Complex_high_prec(Complex_high_prec)> f, Complex_high_prec z0, double tol = 1e-10, int min_inter = 10, int max_iter = 100, double h_eps = 1e-10)
{
    Complex_high_prec z = z0;

    for (int i = 0; i < max_iter; ++i) {
        Complex_high_prec fz = f(z);
        
        //std::cout << i << " " << std::scientific << std::setprecision(15) << "z = " << z << ", f(z) = " << abs(fz) << std::endl;
        
        if (abs(fz) < tol && i > min_inter)
            return z;

        // Approximate derivative f'(z) using complex step
        Complex_high_prec h(h_eps, h_eps);
        Complex_high_prec df = (f(z + h) - fz) / h;

        if (abs(df) < tol)
            throw std::runtime_error("Derivative too small; possible division by zero.");

        z = z - fz / df;
    }
    
    if (abs(z.imag()) < 1e-15 * abs(z.real()))
    {
        std::cout << "Max iterations reached without convergence. The imaginary part is very small, so we go ahead by setting it to zero." << std::endl;
        z.imag(0.0);
        return z;
    }
    
    std::cout << "Max iterations reached without convergence. We go ahead anyway." << std::endl;
    return z;

    //throw std::runtime_error("Max iterations reached without convergence.");
}

struct omegaParams {
    int n;
    int l;
    int m;
    double alpha;
    double atilde;

    // Define less-than operator for map ordering
    bool operator<(const omegaParams& other) const {
        return std::tie(n, l, m, alpha, atilde) <
               std::tie(other.n, other.l, other.m, other.alpha, other.atilde);
    }
};

std::map<omegaParams, std::complex<double>> omega_cache;

// Numerical computation of the complex energy eigenvalue, using Leaver's continued fraction method
std::complex<double> omega_Leaver (int n, int l, int m, double alpha, double atilde)
{
    omegaParams parameters = {n, l, m, alpha, atilde,};
    if (omega_cache.count(parameters))
        return omega_cache[parameters];
    
    int n_max = 2000;
    const Complex_high_prec I(0.0, 1.0);
    double b = sqrt(1-pow(atilde,2));
        
    auto Lambda = [&](Complex_high_prec omega) -> Complex_high_prec
    { return SpheroidalEigenvalue(l, m, pow(atilde,2)*(pow(omega,2) - pow(alpha,2))); };
    auto q = [&](Complex_high_prec omega)
    { return -sqrt(alpha*alpha - omega*omega); };
    
    auto c0 = [&](Complex_high_prec omega) -> Complex_high_prec
        { return 1.0 - 2.0*I*omega - (2.0*I/b) * (omega - atilde*m/2.0); };
    auto c1 = [&](Complex_high_prec omega) -> Complex_high_prec
    { return -4.0 + 4.0*I*(omega - I*q(omega)*(1+b)) + (4.0*I/b) * (omega - atilde*m/2.0) - 2.0*(pow(omega,2) + pow(q(omega),2))/q(omega); };
    auto c2 = [&](Complex_high_prec omega) -> Complex_high_prec
    { return 3.0 - 2.0*I*omega - 2.0*(pow(q(omega),2) - pow(omega,2))/q(omega) - (2.0*I/b) * (omega - atilde*m/2.0); };
    auto c3 = [&](Complex_high_prec omega) -> Complex_high_prec
    { return 2.0*I*pow(omega-I*q(omega),3)/q(omega) + 2.0*pow(omega-I*q(omega),2)*b + pow(q(omega)*atilde,2) + (2.0*m)*I*q(omega)*atilde - Lambda(omega) - 1.0 - pow(omega-I*q(omega),2)/q(omega) + 2.0*q(omega)*b + (2.0*I/b) * (pow(omega-I*q(omega),2)/q(omega) + 1.0) * (omega - atilde*m/2.0); };
    auto c4 = [&](Complex_high_prec omega) -> Complex_high_prec
    { return pow(omega-I*q(omega),4)/pow(q(omega),2) + 2.0*I*omega*pow(omega-I*q(omega),2)/q(omega) - (2.0*I/b) * (pow(omega-I*q(omega),2)/q(omega)) * (omega - atilde*m/2.0); };
    
    auto alpha_n = [&](Complex_high_prec omega, int n) -> Complex_high_prec
    { return pow(n,2) + (c0(omega)+1.0) * (1.0*n) + c0(omega); };
    auto beta_n = [&](Complex_high_prec omega, int n) -> Complex_high_prec
    { return -2.0*pow(n,2) + (c1(omega)+2.0) * (1.0*n) + c3(omega); };
    auto gamma_n = [&](Complex_high_prec omega, int n) -> Complex_high_prec
    { return pow(n,2) + (c2(omega)-3.0) * (1.0*n) + c4(omega); };
    
    auto continued_fraction = [&](Complex_high_prec omega) -> Complex_high_prec
    {
        Complex_high_prec denom = beta_n(omega, n_max);
        for (int n = n_max - 1; n >= 0; --n)
            denom = beta_n(omega, n) - (alpha_n(omega, n) * gamma_n(omega, n+1)) / denom;
        return denom;
    };
    
    Complex_high_prec omega_analytical = alpha + Energy_bound(n, l, m, alpha, atilde) + I * Gamma_Detweiler(n, l, m, alpha, atilde);
    Complex_high_prec omega_high_prec = complex_newton(continued_fraction, omega_analytical);
    
    std::complex<double> omega = std::complex<double>(static_cast<double>(omega_high_prec.real()), static_cast<double>(omega_high_prec.imag()));
    
    omega_cache[parameters] = omega;

    return omega;
}

// Numerical computation of the BH spin on the superradiant threshold
double atilde_threshold (int n, int l, int m, double alpha)
{
    double atilde_analytical = (double) 4 * m * alpha / (pow(m,2) + 4 * pow(alpha,2));
    if (2*alpha > m) atilde_analytical = 0.999;
    
    auto threshold = [&](double atilde) -> double
    { return omega_Leaver(n, l, m, alpha, atilde).real() - m * atilde / (2 * (1 + sqrt(1-pow(atilde,2)))); };
    
    return real_newton(threshold, atilde_analytical);
}

// Kepler's formula, intentionally neglecting the (1+q), as everything else is to leading order in q
double OmegaKepler (double R_star, double q)
{
    return pow(R_star, -3.0/2);
}

struct integralParams {
    int lprime;
    int l_star;
    int l;
    int nprime;
    int n;
    double R_star;
    double alpha;

    // Define less-than operator for map ordering
    bool operator<(const integralParams& other) const {
        return std::tie(lprime, l_star, l, nprime, n, R_star, alpha) <
               std::tie(other.lprime, other.l_star, other.l, other.nprime, other.n, other.R_star, other.alpha);
    }
};

// Inner radial integrand
double innerIntegrand (double r, void * p)
{
    struct integralParams * params = (struct integralParams *)p;
        int lprime = (params->lprime);
        int l_star = (params->l_star);
        int l = (params->l);
        int nprime = (params->nprime);
        int n = (params->n);
        double R_star = (params->R_star);
        double alpha = (params->alpha);
    
    return Rbound(nprime, lprime, r, alpha) * Rbound(n, l, r, alpha) * (0 + pow(r, l_star+2) / pow(R_star, l_star+1));
}

// Outer radial integrand for l_star ≠ 1
double outerIntegrand (double r, void * p)
{
    struct integralParams * params = (struct integralParams *)p;
        int lprime = (params->lprime);
        int l_star = (params->l_star);
        int l = (params->l);
        int nprime = (params->nprime);
        int n = (params->n);
        double R_star = (params->R_star);
        double alpha = (params->alpha);
    
    return Rbound(nprime, lprime, r, alpha) * Rbound(n, l, r, alpha) * (0 + pow(R_star, l_star) / pow(r, l_star-1));
}

// Outer radial integrand for l_star = 1
double outerIntegrand_dipole (double r, void * p)
{
    struct integralParams * params = (struct integralParams *)p;
        int lprime = (params->lprime);
        int l_star = (params->l_star);
        int l = (params->l);
        int nprime = (params->nprime);
        int n = (params->n);
        double R_star = (params->R_star);
        double alpha = (params->alpha);
    
    return Rbound(nprime, lprime, r, alpha) * Rbound(n, l, r, alpha) * (- pow(r,l_star)/pow(R_star,l_star+1) + pow(R_star,l_star)/pow(r,l_star+1)) * pow(r,2);
}

std::map<integralParams, double> I_r_cache;

// Radial integral for l_star ≠ 1
double I_r (gsl_integration_workspace * w, int lprime, int l_star, int l, int nprime, int n, double R_star, double alpha)
{
    assert (l_star != 1);
    
    int max_subdivisions = 1000;
    
    double relativeError = 1e-6;
        
    double innerIntegral, innerError, outerIntegral, outerError;
    integralParams parameters = {lprime, l_star, l, nprime, n, R_star, alpha,};
    if (I_r_cache.count(parameters))
        return I_r_cache[parameters];
    
    gsl_function Inner;
    Inner.function = &innerIntegrand;
    Inner.params = &parameters;
    
    gsl_integration_qag(&Inner, 0, R_star, 0, relativeError, max_subdivisions, 6, w, &innerIntegral, &innerError);
    
    gsl_function Outer;
    Outer.function = &outerIntegrand;
    Outer.params = &parameters;
    
    gsl_integration_qagiu(&Outer, R_star, 0, relativeError, max_subdivisions, w, &outerIntegral, &outerError);
    
    // When R_star is too large, outerIntegral = nan, which happens because there is a huge exponential suppression, so we just manually set it to zero.
    if(innerIntegral!=innerIntegral) {/*std::cout << "warning! 1" << std::endl;*/ innerIntegral = 0;}
    if(outerIntegral!=outerIntegral) {/*std::cout << "warning! 2" << std::endl;*/ outerIntegral = 0;}
        
    I_r_cache[parameters] = innerIntegral + outerIntegral;
    return innerIntegral + outerIntegral;
}

// Radial integral for l_star = 1
double I_r_dipole (gsl_integration_workspace * w, int lprime, int l_star, int l, int nprime, int n, double R_star, double alpha)
{
    assert (l_star == 1);
    
    int max_subdivisions = 1000;
    
    double relativeError = 1e-6;
        
    double outerIntegral, outerError;
    integralParams parameters = {lprime, l_star, l, nprime, n, R_star, alpha,};
    if (I_r_cache.count(parameters))
        return I_r_cache[parameters];
    
    gsl_function Outer;
    Outer.function = &outerIntegrand_dipole;
    Outer.params = &parameters;
    
    gsl_integration_qagiu(&Outer, R_star, 0, relativeError, max_subdivisions, w, &outerIntegral, &outerError);
    
    // When R_star is too large, outerIntegral = nan, which happens because there is a huge exponential suppression, so we just manually set it to zero.
    if(outerIntegral!=outerIntegral) {/*std::cout << "warning! 3" << std::endl;*/ outerIntegral = 0;}
    
    I_r_cache[parameters] = outerIntegral;
    return outerIntegral;
}

// Matrix element of V_star, as a functoon of the eccentric anomaly E
std::complex<double> V_star_matrix_element (gsl_integration_workspace * w, int nprime, int lprime, int mprime, int n, int l, int m, double a, int g, double alpha, double q, double eccentricity, double inclination, double argument_periapsis, double E)
{
    const std::complex<double> I(0.0, 1.0);
    
    double true_anomaly;
    
    if (E < M_PI) true_anomaly = 2 * atan( sqrt( (1+eccentricity)*pow(tan(E/2),2)/(1-eccentricity) ) );
    else true_anomaly = 2 * M_PI - 2 * atan( sqrt( (1+eccentricity)*pow(tan(E/2),2)/(1-eccentricity) ) );
    
    double R_star = a * (1 - eccentricity * cos(E));
    
    double radial_integral;
    std::complex<double> sph_harmonic = 0, result = 0;
    
    for (int l_star = std::max(abs(m-mprime),abs(lprime-l)); l_star <= l+lprime; l_star++)
    {
        if (l_star != 1)
            radial_integral = I_r(w, lprime, l_star, l, nprime, n, R_star, alpha);
        else
            radial_integral = I_r_dipole(w, lprime, l_star, l, nprime, n, R_star, alpha);
        
        sph_harmonic = 0;
        for (int m_tilde = -l_star; m_tilde <= l_star; m_tilde++)
            sph_harmonic += Wigner_small_d(l_star, m-mprime, -m_tilde, inclination) * gsl_sf_legendre_sphPlm(l_star, abs(-m_tilde), 0) * std::exp(- I * (double)m_tilde * (true_anomaly - argument_periapsis)) / (2.0*l_star + 1);
                
        result += sph_harmonic * radial_integral * I_Omega(lprime, l_star, l, mprime, m-mprime, m);
    }
        
    return result * (- 4 * M_PI * alpha * q);
}

// Trapezoidal integration for complex functions. Note that N must be much larger than g to avoid undersampling
std::complex<double> trapezoidal_integration(const std::function<std::complex<double>(double)>& f, double a, double b, int N = 100)
{
    double h = (b - a) / N;
    std::complex<double> sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < N; ++i) {
        sum += f(a + i * h);
    }
    return h * sum;
}

struct eta_g_cache_struct { std::vector<double> a_grid; std::vector<double> values; };

using Key = std::tuple<int, int, int, int, int, int, int, double>; // (n', l', m', n, l, m, g, alpha)

std::map<Key, eta_g_cache_struct> eta_g_cache;

// Complex Fourier coefficient of the matrix element of the perturbation
double eta_g(gsl_integration_workspace * w, int nprime, int lprime, int mprime, int n, int l, int m, double a, int g, double alpha, double q, double eccentricity, double inclination, double argument_periapsis)
{
    const std::complex<double> I(0.0, 1.0);
    auto integrand = [&](double E) {
        return V_star_matrix_element(w, nprime, lprime, mprime, n, l, m, a, g, alpha, q, eccentricity, inclination, argument_periapsis, E) * std::exp(I * (double)g * (E - eccentricity*sin(E))) * (1 - eccentricity*cos(E));
    };
        
    return abs(trapezoidal_integration(integrand, 0.0, 2*M_PI) / (2*M_PI));
}

// Complex Fourier coefficient of the matrix element of the perturbation (optimized for equatorial circular orbits)
double eta_g_circular_equatorial (gsl_integration_workspace * w, int nprime, int lprime, int mprime, int n, int l, int m, double R_star, double alpha, double q)
{
    double radial_integral;
    double sph_harmonic = 0, result = 0;
    
    for (int l_star = std::max(abs(m-mprime),abs(lprime-l)); l_star <= l+lprime; l_star++)
    {
        if (l_star != 1)
            radial_integral = I_r(w, lprime, l_star, l, nprime, n, R_star, alpha);
        else
            radial_integral = I_r_dipole(w, lprime, l_star, l, nprime, n, R_star, alpha);
        
        sph_harmonic = gsl_sf_legendre_sphPlm(l_star, abs(m-mprime), 0) / (2.0*l_star + 1);
                
        result += sph_harmonic * radial_integral * I_Omega(lprime, l_star, l, mprime, m-mprime, m);
    }
                        
    return result * (- 4 * M_PI * alpha * q);
}

// Pre-compute and store the values of eta_g. These values will then be linearly interpolated when computing the horizon fluxes. This is a very significant optimization when the horizon fluxes need to be computed on a fine grid in order to resolve the narrow resonant peaks
void precompute_eta_g_cache (gsl_integration_workspace* w, int n, int l, int m, const std::vector<double> & a_grid, int nprime_max, int lprime_max, int g_max, double alpha, double q, double eccentricity, double inclination, double argument_periapsis)
{
    for (int nprime = 1; nprime <= nprime_max; ++nprime)
        for (int lprime = 0; lprime <= std::min(nprime - 1, lprime_max); ++lprime)
            for (int mprime = -lprime; mprime <= lprime; ++mprime)
            {
                for (int g = -g_max; g <= g_max; g++)
                {
                    if (nprime == n && lprime == l && mprime == m) continue;
                    std::cout << "Computing eta on a grid for: " << nprime << " " << lprime << " " << mprime << ", g = " << g << std::endl;
                    Key key = {nprime, lprime, mprime, n, l, m, g, alpha};
                    eta_g_cache_struct cache;
                    cache.a_grid = a_grid;
                    for (double a : a_grid)
                        cache.values.push_back(eta_g(w, nprime, lprime, mprime, n, l, m, a, g, alpha, q, eccentricity, inclination, argument_periapsis));
                    eta_g_cache[key] = std::move(cache);
                }
            }
}

// Pre-compute and store the values of eta_g (optimized for equatorial circular orbits). These values will then be linearly interpolated when computing the horizon fluxes. This is a very significant optimization when the horizon fluxes need to be computed on a fine grid in order to resolve the narrow resonant peaks
void precompute_eta_g_circular_equatorial_cache (gsl_integration_workspace* w, int n, int l, int m, int sign_orbit, const std::vector<double> & R_grid, int nprime_max, int lprime_max, double alpha, double q)
{
    for (int nprime = 1; nprime <= nprime_max; ++nprime)
        for (int lprime = 0; lprime <= std::min(nprime - 1, lprime_max); ++lprime)
            for (int mprime = -lprime; mprime <= lprime; ++mprime)
                
            {
                std::cout << "Computing eta on a grid for: " << nprime << " " << lprime << " " << mprime << std::endl;
                if (nprime == n && lprime == l && mprime == m) continue;
                Key key = {nprime, lprime, mprime, n, l, m, sign_orbit*(mprime-m), alpha};
                eta_g_cache_struct cache;
                cache.a_grid = R_grid;
                for (double R : R_grid)
                    cache.values.push_back(eta_g_circular_equatorial(w, nprime, lprime, mprime, n, l, m, R, alpha, q));
                eta_g_cache[key] = std::move(cache);
            }
}

// Just a linear interpolator
double linear_interp(const std::vector<double>& x, const std::vector<double>& y, double x_new)
{
    assert(x.size() == y.size());
    auto it = std::lower_bound(x.begin(), x.end(), x_new);

    if (it == x.begin()) return y.front();
    if (it == x.end()) return y.back();

    auto idx = std::distance(x.begin(), it);
    double x0 = x[idx - 1], x1 = x[idx];
    double y0 = y[idx - 1], y1 = y[idx];

    double t = (x_new - x0) / (x1 - x0);
    return y0 * (1 - t) + y1 * t;
}

struct HorizonFluxes {std::vector<double> R_grid; std::vector<double> HorizonMassRate; std::vector<double> HorizonPower; std::vector<double> HorizonTorque_z; std::vector<double> HorizonTorque_x;};

// Create a grid for the semi-major axis that is designed to properly sample the (potentially extremely narrow) resonance peak
std::vector<double> R_grid_resonance (double DeltaE, double DeltaGamma, int g, double Rmin, double Rmax)
{
    const int N_1 = 101;
    std::vector<double> tan_values;
    double Omega, Omegamin = pow(Rmax,-1.5), Omegamax = pow(Rmin,-1.5);
    if (abs(DeltaGamma) < 1e-35) DeltaGamma = -1e-35; // Safety fallback in case DeltaGamma = 0

    for (int i = 1; i <= N_1; i++)
    {
        Omega = (std::tan(- M_PI / 2 + i * M_PI / (N_1 + 1)) * abs(DeltaGamma) + DeltaE) / g ;
        if (Omega > Omegamin && Omega < Omegamax) tan_values.push_back( pow(Omega, -2.0/3) );
        if (Omega > Omegamin && Omega < Omegamax && i == 1 + N_1/2) std::cout << "Possible resonance with g = " << g << ", at R_res = " << pow(Omega, -2.0/3) << std::endl;
    }
        
    const int N_2 = Rmax-Rmin;
    for (double R = Rmin; R <= Rmax; R += (Rmax-Rmin)/N_2)
        tan_values.push_back(R);
    
    double epsilon = 0.1;
    int N = log( pow((Rmax-Rmin)/N_2, -1.5) / (abs(DeltaGamma) * epsilon) ) / epsilon;
    for (int i = - N/3; i <= N; i++)
    {
        Omega = (std::exp(epsilon * i) * abs(DeltaGamma) + DeltaE) / g ;
        if (Omega > Omegamin && Omega < Omegamax) tan_values.push_back( pow(Omega, -2.0/3) );
        Omega = (- std::exp(epsilon * i) * abs(DeltaGamma) + DeltaE) / g ;
        if (Omega > Omegamin && Omega < Omegamax) tan_values.push_back( pow(Omega, -2.0/3) );
    }
        
    std::sort(tan_values.begin(), tan_values.end());
        
    return tan_values;
}

// Compute the horizon fluxes
HorizonFluxes horizon_flux(gsl_integration_workspace * w, int n, int l, int m, int nprime_max, int lprime_max, int g_max, double alpha, double McOverM, double q, double eccentricity, double inclination, double argument_periapsis, double Rmin, double Rmax)
{
    double atilde;
    if (m == 0)
        atilde = 0;
    else
        atilde = atilde_threshold(n, l, m, alpha);
    std::complex<double> omega_a = omega_Leaver(n, l, m, alpha, atilde);
    double E_a = omega_a.real();
    double Gamma_a = 0 * omega_a.real(); // manually setting it to zero because numerical errors on Gamma_a can be larger than |Gamma_b|
    
    std::complex<double> omega_b;
    double E_b, Gamma_b, tmp;
            
    std::vector<HorizonFluxes> List_resonance_fluxes;
    
    for (int nprime = 1; nprime <= nprime_max; nprime++)
        for (int lprime = 0; lprime <= std::min(nprime-1, lprime_max); lprime++)
            for (int mprime = -lprime; mprime <= lprime; mprime++)
                {
                    if (nprime == n && lprime == l && mprime == m) continue;
                    //if (!(lprime == 2 && mprime == -2)) continue; // uncomment if want to study only one mode
                    
                    omega_b = omega_Leaver(nprime, lprime, mprime, alpha, atilde);
                    E_b = omega_b.real();
                    Gamma_b = omega_b.imag();
                    
                    for (int g = -g_max; g <= g_max; g++)
                    {
                        std::cout << "Computing the flux for " << nprime << " " << lprime << " " << mprime << ", g = " << g << std::endl;
                        HorizonFluxes fluxes;
                        Key key = {nprime, lprime, mprime, n, l, m, g, alpha};
                        
                        fluxes.R_grid = R_grid_resonance(E_b-E_a, Gamma_b-Gamma_a, g, Rmin, Rmax);
                        
                        for (double R_star : fluxes.R_grid)
                        {
                            double Omega = OmegaKepler(R_star, q);
                            double eta_val = linear_interp(eta_g_cache.at(key).a_grid, eta_g_cache.at(key).values, R_star);
                            tmp = 2 * (Gamma_b - Gamma_a) * pow(abs(eta_val), 2) / (pow(E_a - E_b + g*Omega, 2) + pow(Gamma_b - Gamma_a, 2));
                            
                            fluxes.HorizonMassRate.push_back(tmp);
                            fluxes.HorizonPower.push_back(- McOverM * (E_b - E_a) * tmp / alpha);
                            fluxes.HorizonTorque_z.push_back(- McOverM * (mprime - m) * tmp / alpha);
                        }
                        
                        List_resonance_fluxes.push_back(fluxes);
                    }
                }
    
    
    HorizonFluxes result;
    
    std::set<double> R_union_set;
    for (const auto& fluxes : List_resonance_fluxes)
        R_union_set.insert(fluxes.R_grid.begin(), fluxes.R_grid.end());
    
    result.R_grid.assign(R_union_set.begin(), R_union_set.end());
    
    result.HorizonMassRate.assign(result.R_grid.size(), 0.0);
    result.HorizonPower.assign(result.R_grid.size(), 0.0);
    result.HorizonTorque_z.assign(result.R_grid.size(), 0.0);
        
    for (const auto & fluxes : List_resonance_fluxes)
    {
        for (size_t i = 0; i < result.R_grid.size(); ++i)
        {
            double R_star = result.R_grid[i];
            result.HorizonMassRate[i] += linear_interp(fluxes.R_grid, fluxes.HorizonMassRate, R_star);
            result.HorizonPower[i] += linear_interp(fluxes.R_grid, fluxes.HorizonPower, R_star);
            result.HorizonTorque_z[i] += linear_interp(fluxes.R_grid, fluxes.HorizonTorque_z, R_star);
        }
    }
        
    return result;
}

// Compute the horizon fluxeso (ptimized for equatorial circular orbits)
HorizonFluxes horizon_flux_circular_equatorial(gsl_integration_workspace * w, int n, int l, int m, int sign_orbit, int nprime_max, int lprime_max, double alpha, double McOverM, double q, double Rmin, double Rmax)
{
    double atilde;
    if (m == 0)
        atilde = 0;
    else
        atilde = atilde_threshold(n, l, m, alpha);
    std::complex<double> omega_a = omega_Leaver(n, l, m, alpha, atilde);
    double E_a = omega_a.real();
    double Gamma_a = 0 * omega_a.real(); // setting it to zero because numerical errors can be larger than Gamma_b
    
    std::complex<double> omega_b;
    double E_b, Gamma_b, tmp;
                
    std::vector<HorizonFluxes> List_resonance_fluxes;
    
    for (int nprime = 1; nprime <= nprime_max; nprime++)
        for (int lprime = 0; lprime <= std::min(nprime-1, lprime_max); lprime++)
            for (int mprime = -lprime; mprime <= lprime; mprime++)
                {
                    if (nprime == n && lprime == l && mprime == m) continue;
                    //if (!(lprime == 2 && mprime == -2)) continue; // uncomment if want to study only one mode
                    
                    std::cout << "Computing the flux for " << nprime << " " << lprime << " " << mprime << std::endl;

                    HorizonFluxes fluxes;
                    omega_b = omega_Leaver(nprime, lprime, mprime, alpha, atilde);
                    E_b = omega_b.real();
                    Gamma_b = omega_b.imag();
                    Key key = {nprime, lprime, mprime, n, l, m, sign_orbit*(mprime-m), alpha};
                    
                    fluxes.R_grid = R_grid_resonance(E_b-E_a, Gamma_b-Gamma_a, (mprime-m)*sign_orbit, Rmin, Rmax);
                    
                    for (double R_star : fluxes.R_grid)
                    {
                        double Omega = OmegaKepler(R_star, q);
                        double eta_val = linear_interp(eta_g_cache.at(key).a_grid, eta_g_cache.at(key).values, R_star);
                        tmp = 2 * (Gamma_b - Gamma_a) * pow(abs(eta_val), 2) / (pow(E_a - E_b + sign_orbit*(mprime-m)*Omega, 2) + pow(Gamma_b - Gamma_a, 2));
                        
                        fluxes.HorizonMassRate.push_back(tmp);
                        fluxes.HorizonPower.push_back(- McOverM * (E_b - E_a) * tmp / alpha);
                        fluxes.HorizonTorque_z.push_back(- McOverM * (mprime - m) * tmp / alpha);
                    }
                    
                    List_resonance_fluxes.push_back(fluxes);
                }
    
    
    HorizonFluxes result;
    
    std::set<double> R_union_set;
    for (const auto& fluxes : List_resonance_fluxes)
        R_union_set.insert(fluxes.R_grid.begin(), fluxes.R_grid.end());
    
    result.R_grid.assign(R_union_set.begin(), R_union_set.end());
    
    result.HorizonMassRate.assign(result.R_grid.size(), 0.0);
    result.HorizonPower.assign(result.R_grid.size(), 0.0);
    result.HorizonTorque_z.assign(result.R_grid.size(), 0.0);
        
    for (const auto & fluxes : List_resonance_fluxes)
    {
        for (size_t i = 0; i < result.R_grid.size(); ++i)
        {
            double R_star = result.R_grid[i];
            result.HorizonMassRate[i] += linear_interp(fluxes.R_grid, fluxes.HorizonMassRate, R_star);
            result.HorizonPower[i] += linear_interp(fluxes.R_grid, fluxes.HorizonPower, R_star);
            result.HorizonTorque_z[i] += linear_interp(fluxes.R_grid, fluxes.HorizonTorque_z, R_star);
        }
    }
        
    return result;
}

int main () {
    
    gsl_set_error_handler_off();
        
    int max_subdivisions_K = 1000;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (max_subdivisions_K);
    
    int nprime_max = 2, lprime_max = 1, abs_g_max = 20;
    int n = 2, l = 1, m = 1;
    double alpha = 0.2;
    double McOverM = 0.01;
    double q = 0.001;
    double eccentricity = 0.5, inclination = 70*M_PI/180, argument_periapsis = 0*M_PI/180;
    
    HorizonFluxes Flux;
    
    // Range of semi-major axes
    double Rmin = 1000, Rmax = 80000;
    
    std::vector<double> R_grid;
    for (double R = Rmin; R <= Rmax; R += R*0.01) R_grid.push_back(R);
    
    int sign_orbit = -1; // 1 for prograde, -1 for retrograde
    
    // -------------------------------------------------------------------
    // Uncomment the following lines to compute the horizon fluxes on an equatorial orbit:
    /*
    precompute_eta_g_circular_equatorial_cache(w, n, l, m, sign_orbit, R_grid, nprime_max, lprime_max, alpha, q);
        
    Flux = horizon_flux_circular_equatorial(w, n, l, m, sign_orbit, nprime_max, lprime_max, alpha, McOverM, q, Rmin, Rmax);
    for (size_t i = 0; i < Flux.R_grid.size(); ++i)
        std::cout << std::setprecision(25) << Flux.R_grid[i] << " " << Flux.HorizonMassRate[i] << " " << Flux.HorizonPower[i] << std::endl;
    */
    // -------------------------------------------------------------------
    
    // -------------------------------------------------------------------
    // Uncomment the following lines to compute the horizon fluxes on an orbit with eccentricity, inclination and argument_periapsis specified above:
    
    precompute_eta_g_cache(w, n, l, m, R_grid, nprime_max, lprime_max, abs_g_max, alpha, q, eccentricity, inclination, argument_periapsis);
    Flux = horizon_flux(w, n, l, m, nprime_max, lprime_max, abs_g_max, alpha, McOverM, q, eccentricity, inclination, argument_periapsis, Rmin, Rmax);
    for (size_t i = 0; i < Flux.R_grid.size(); ++i)
        std::cout << std::setprecision(25) << Flux.R_grid[i] << " " << Flux.HorizonMassRate[i] << " " << Flux.HorizonPower[i] << std::endl;
    // -------------------------------------------------------------------
    
    // -------------------------------------------------------------------
    // Uncomment the following lines to compute Pres for various values of alpha:
    /*
    std::cout << std::setprecision(25);
    for (double alpha = 0.05; alpha <= 0.35; alpha += 0.001)
    {
        precompute_eta_g_cache(w, n, l, m, R_grid, nprime_max, lprime_max, g_max, alpha, q, eccentricity, inclination, argument_periapsis);
        Flux = horizon_flux(w, n, l, m, nprime_max, lprime_max, g_max, alpha, McOverM, q, eccentricity, inclination, argument_periapsis, Rmin, Rmax);
        std::cout << alpha << " " << Flux.R_grid.size() << " ";
        for (size_t i = 0; i < Flux.R_grid.size(); ++i)
            std::cout << Flux.R_grid[i] << " ";
        for (size_t i = 0; i < Flux.R_grid.size(); ++i)
            std::cout << Flux.HorizonPower[i] << " ";
        std::cout << std::endl;
    }
    */
    // -------------------------------------------------------------------
    
    gsl_integration_workspace_free(w);
    
    return 0;
}
