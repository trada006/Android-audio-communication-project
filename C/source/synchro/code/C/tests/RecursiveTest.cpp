/* 
 * File:   RecursiveTest.cpp
 * Author: Robby McKilliam
 *
 * Created on 26/01/2013, 2:03:19 PM
 */

#include <cstdlib>
#include <iostream>
#include <complex>
#include "FastAlgorithm.h"
#include "TabulatedFastAlgorithm.h"
#include "FinitePulse.h"

using namespace std;

static double tol = 1e-3;
static double T = 1.0; //symbol period
static unsigned int p = 1;
static unsigned int q = 5;
static unsigned int c = 4; //new c, its will divide T
static double Ts = p * T / q; //sample period
static double Delta = T / c;
static double taumin = 10.0;
static double taumax = 30.0;
static double tau0 = 20.111;
static unsigned int L = 100;
static unsigned int N = (int) ceil((T * L + taumax) / Ts);
static int numzeros = 2;
static vector<int> P;
static vector<int> D;
static vector<complex<double> > pilots;
static vector<complex<double> > s;
static vector<complex<double> > r;
static unsigned int M = 4; //QPSK
static int numpilots = 40;
static TruncatedSincPulse tx(T, numzeros);
static complex<double> a0 = polar<double>(1,0.1*pi); //phase and amplitude

void testRecursive() {

    Recursive<TruncatedSincPulse> rec(P, D, pilots, tx, T, Ts, taumin, taumax, c, p, q);
    Naive<TruncatedSincPulse> nai(P, D, pilots, tx, T, Ts, taumin, taumax, c);
    double tauhatrec = rec.estimate(r);
    double tauhatnai = nai.estimate(r);

    std::cout << "Recursive Zgrid test ... ";
    bool pass = true;
    vector<double> zgrid = rec.getZgrid();
    for (int i = 0; i < zgrid.size(); i++) {
        pass &= fabs(zgrid[i] - nai.Z(taumin + i * Delta)) < tol;
        pass &= fabs(zgrid[i] - rec.Z(taumin + i * Delta)) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;

    std::cout << "Recursive Ygrid test ... ";
    pass = true;
    vector<double> ygrid = rec.getYgrid();
    for (int i = 0; i < ygrid.size(); i++) {
        pass &= fabs(ygrid[i] - std::abs(nai.Y(taumin + i * Delta))) < tol;
        pass &= fabs(ygrid[i] - std::abs(rec.Y(taumin + i * Delta))) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;

    std::cout << "Recursive SSgrid test ... ";
    pass = true;
    vector<double> ssgrid = rec.getSSgrid();
    for (int i = 0; i < ssgrid.size(); i++) {
        pass &= fabs(ssgrid[i] - std::abs(nai.SS(taumin + i * Delta))) < tol;
        pass &= fabs(ssgrid[i] - std::abs(rec.SS(taumin + i * Delta))) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;

    std::cout << "Recursive coarse estimate... ";
    pass = fabs(rec.coarseMaximiseSS() - nai.coarseMaximiseSS()) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "Recursive fine estimate... ";
    pass = fabs(tauhatrec - tauhatnai) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;

}

void testTabulatedRecursive() {

    TabulatedRecursive<TruncatedSincPulse> rec(P, D, pilots, tx, T, Ts, taumin, taumax, c, p, q);
    Naive<TruncatedSincPulse> nai(P, D, pilots, tx, T, Ts, taumin, taumax, c);
    double tauhatrec = rec.estimate(r);
    double tauhatnai = nai.estimate(r);

    std::cout << "TabulatedRecursive Zgrid test ... ";
    bool pass = true;
    vector<double> zgrid = rec.getZgrid();
    for (int i = 0; i < zgrid.size(); i++) {
        pass &= fabs(zgrid[i] - nai.Z(taumin + i * Delta)) < tol;
        pass &= fabs(zgrid[i] - rec.Z(taumin + i * Delta)) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;

    std::cout << "TabulatedRecursive Ygrid test ... ";
    pass = true;
    vector<double> ygrid = rec.getYgrid();
    for (int i = 0; i < ygrid.size(); i++) {
        pass &= fabs(ygrid[i] - std::abs(nai.Y(taumin + i * Delta))) < tol;
        pass &= fabs(ygrid[i] - std::abs(rec.Y(taumin + i * Delta))) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;

    std::cout << "TabulatedRecursive SSgrid test ... ";
    pass = true;
    vector<double> ssgrid = rec.getSSgrid();
    for (int i = 0; i < ssgrid.size(); i++) {
        pass &= fabs(ssgrid[i] - std::abs(nai.SS(taumin + i * Delta))) < tol;
        pass &= fabs(ssgrid[i] - std::abs(rec.SS(taumin + i * Delta))) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;

    std::cout << "TabulatedRecursive coarse estimate... ";
    pass = fabs(rec.coarseMaximiseSS() - nai.coarseMaximiseSS()) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "TabulatedRecursive fine estimate... ";
    pass = fabs(tauhatrec - tauhatnai) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;

}

int main(int argc, char** argv) {
    
    for (int m = 0; m < numpilots; m++) P.push_back(m);
    for (int m = numpilots; m < L; m++) D.push_back(m);
    for (int m = 1; m <= L; m++) s.push_back(polar<double>(1.0, (2*pi*(rand()%M))/M));
    for (int i = 0; i < P.size(); i++) pilots.push_back(s[P[i]]);
    
    //function generates transmitted signal
    auto x = [] (double t){
        int mini = max((int) 1, (int) ceil((t - tx.tmax()) / T));
        int maxi = min((int) L, (int) floor((t - tx.tmin()) / T));
        complex<double> sum(0, 0);
        for (int i = mini; i <= maxi; i++) sum += s[i - 1] * tx.pulse(t - i * T); //this here is already fairly time consuming!
        return sum;
    };
     for(int n = 1; n <= N; n++) r.push_back(a0*x(n*Ts-tau0)); //generate noiseless received signal
    
    
    testRecursive();
    testTabulatedRecursive();

    return 0;
}
