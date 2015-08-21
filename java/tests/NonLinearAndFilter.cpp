/* 
 * File:   NonLinearAndFilter.cpp
 * Author: mckillrg
 *
 * Created on 12/03/2013, 10:17:17 AM
 */

#include <cstdlib>
#include <iostream>
#include <complex>
#include "SymbolTiming.h"
#include "FinitePulse.h"

using namespace std;

static double tol = 1e-2;
static double T = 1.0; //symbol period
static unsigned int p = 3;
static unsigned int q = 21;
static unsigned int c = 11; //new c, its will divide T
static double Ts = p * T / q; //sample period
static double Delta = T/c;
static double taumin = 10.0;
static double taumax = 30.0;
static double tau0 = 20.111;
static unsigned int L = 100;
static unsigned int N = (int) ceil((T * L + taumax) / Ts);
static int numzeros = 2;
static vector<int> P;
static vector<int> D;
static vector<int> S;
static vector<complex<double> > pilots;
static vector<complex<double> > s;
static vector<complex<double> > r;
static unsigned int M = 4; //QPSK
static int numpilots = 30;
static TruncatedRootRaisedCosine tx(T, 1.0/10.0, numzeros);
static complex<double> a0 = polar<double>(1,0.1*pi); //phase and amplitude

void testSquareLaw() {

    auto f = [] (complex<double> x) { return std::norm(x); }; //square law estimator
    NonlinearityAndFilter<TruncatedRootRaisedCosine> est(S, tx, T, Ts, taumin, taumax, f);
    double gammahat = est.estimate(r);
    double gamma0 = tau0 - T*round(tau0/T);
    
    std::cout << "Square law test ... ";
    //std::cout << gammahat << ", " << gamma0 << std::endl;
    bool pass = fabs(gammahat - gamma0) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    TabulatedNonlinearityAndFilter<TruncatedRootRaisedCosine> tabest(S, tx, T, Ts, taumin, taumax, f);
    double tabgammahat = tabest.estimate(r);
    std::cout << "Tabulated Square law test ... ";
    //std::cout << gammahat << ", " << tabgammahat << std::endl;
    pass = fabs(gammahat - tabgammahat) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
}

void testAbsLaw() {
    auto f = [] (complex<double> x) { return std::abs(x); }; //square law estimator
    NonlinearityAndFilter<TruncatedRootRaisedCosine> est(S, tx, T, Ts, taumin, taumax, f);
    double gammahat = est.estimate(r);
    double gamma0 = tau0 - T*round(tau0/T);
    
    std::cout << "Absolute value law test ... ";
    //std::cout << gammahat << ", " << gamma0 << std::endl;
    bool pass = fabs(gammahat - gamma0) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    TabulatedNonlinearityAndFilter<TruncatedRootRaisedCosine> tabest(S, tx, T, Ts, taumin, taumax, f);
    double tabgammahat = tabest.estimate(r);
    std::cout << "Tabulated Absolute law test ... ";
    //std::cout << gammahat << ", " << tabgammahat << std::endl;
    pass = fabs(gammahat - tabgammahat) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
}

int main(int argc, char** argv) {
    
    for (int m = 0; m < numpilots; m++) P.push_back(m);
    for (int m = numpilots; m < L; m++) D.push_back(m);
    for (int m = 0; m < L; m++) S.push_back(m); //contatenate P and D
    for (int m = 1; m <= L; m++) s.push_back(polar<double>(1.0, (2*pi*(rand()%M))/M));
    for ( int i = 0; i < P.size(); i++ ) pilots.push_back(s[P[i]]);

    //function generates transmitted signal
    auto x = [&s, &tx, T, L] (double t){
        int mini = max((int) 1, (int) ceil((t - tx.tmax()) / T));
        int maxi = min((int) L, (int) floor((t - tx.tmin()) / T));
        complex<double> sum(0, 0);
        for (int i = mini; i <= maxi; i++) sum += s[i - 1] * tx.pulse(t - i * T); //this here is already fairly time consuming!
        return sum;
    };
    for(int n = 1; n <= N; n++) r.push_back(a0*x(n*Ts-tau0)); //generate noiseless received signal

    testSquareLaw();
    testAbsLaw();
    
    return 0;
}

