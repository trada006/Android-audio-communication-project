/* 
 * File:   QAMTest.cpp
 * Author: Robby McKilliam
 *
 * Created on 27 January 2015, 10:20 AM
 */

#include <cstdlib>
#include <iostream>
#include <complex>
#include "TimeOffsetEstimator.h"
#include "QAMTimeOffsetEstimator.h"
#include "FinitePulse.h"

using namespace std;

static double tol = 1e-3;
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
static vector<complex<double> > pilots;
static vector<complex<double> > s;
static vector<complex<double> > r;
static unsigned int M = 4; //QPSK
static int numpilots = 30;
static TruncatedSincPulse tx(T, numzeros);
static complex<double> a0 = polar<double>(1,0.1*pi); //phase and amplitude

void testQAMDirect() {

    QAMDirect<TruncatedSincPulse> dir(P, D, pilots, tx, T, Ts, taumin, taumax, c, p, q);
    Naive<TruncatedSincPulse> nai(P, D, pilots, tx, T, Ts, taumin, taumax, c);
    double tauhatdir = dir.estimate(r);
    double tauhatnai = nai.estimate(r);
    
    //These are QPSK symbls so the energy in the pilot sequence should be equal to the number of pilots
    if(std::abs(P.size() - dir.Ep) < 1e-8) std::cout << "QAM Direct Ep ... PASS" << std::endl;
    else std::cout << "QAM Direct Ep ... FAIL" << std::endl;
        
    std::cout << "QAMDirect Zgrid test ... ";
    bool pass = true;
    vector<double> zgrid = dir.getZgrid();
    for(int i = 0; i < zgrid.size(); i++ ) {
        //std::cout << zgrid[i] << ", " << nai.Z(taumin + i*Delta) << std::endl;
        pass &= fabs(zgrid[i] - dir.Z(taumin + i*Delta)) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "QAMDirect Ygrid test ... ";
    pass = true;
    vector<double> ygrid = dir.getYgrid();
    for(int i = 0; i < ygrid.size(); i++ ) {
        pass &= fabs(ygrid[i] - std::norm(dir.Y(taumin + i*Delta))/dir.Ep) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "QAMDirect SSgrid test ... ";
    pass = true;
    vector<double> ssgrid = dir.getSSgrid();
    for(int i = 0; i < ssgrid.size(); i++ ) {
        //std::cout << ssgrid[i] << ", " << std::abs(nai.SS(taumin + i*Delta)) << ", " << std::abs(dir.SS(taumin + i*Delta)) << std::endl;
        pass &= fabs(ssgrid[i] - std::abs(dir.SS(taumin + i*Delta))) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "QAMDirect coarse estimate... ";
    //std::cout << dir.coarseMaximiseSS() << ", " << nai.coarseMaximiseSS() << std::endl;
    pass = fabs(dir.coarseMaximiseSS() - nai.coarseMaximiseSS()) < tol; //these are not the same estimator, but should be close in this noiseless test
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "QAMDirect fine estimate... ";
    //std::cout << tauhatdir << ", " << tauhatnai << std::endl;
    pass = fabs(tauhatdir - tauhatnai) < 1e-1; //these are not the same estimator, but should be close in this noiseless test
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
}

int main(int argc, char** argv) {
    
    for (int m = 0; m < numpilots; m++) P.push_back(m);
    for (int m = numpilots; m < L; m++) D.push_back(m);
    for (int m = 1; m <= L; m++) s.push_back(polar<double>(1.0, (2*pi*(rand()%M))/M));
    for ( int i = 0; i < P.size(); i++ ) pilots.push_back(s[P[i]]);

    //function generates transmitted signal
    auto x = [] (double t){
        int mini = max((int) 1, (int) ceil((t - tx.tmax()) / T));
        int maxi = min((int) L, (int) floor((t - tx.tmin()) / T));
        complex<double> sum(0, 0);
        for (int i = mini; i <= maxi; i++) sum += s[i - 1] * tx.pulse(t - i * T); //this here is already fairly time consuming!
        return sum;
    };
    for(int n = 1; n <= N; n++) r.push_back(a0*x(n*Ts-tau0)); //generate noiseless received signal

    testQAMDirect();
    
    return 0;
}


