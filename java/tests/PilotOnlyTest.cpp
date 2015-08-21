/* 
 * File:   PilotOnlyTest.cpp
 * Author: harprobey
 *
 * Created on 05/03/2013, 10:12:40 PM
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
static unsigned int p = 3;
static unsigned int q = 21;
static unsigned int c = 11; //new c, its will divide T
static double Ts = p * T / q; //sample period
static double Delta = T/c;
static double taumin = 10.0;
static double taumax = 30.0;
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

void test1() {
    Direct<TruncatedSincPulse> dir(P, D, pilots, TruncatedSincPulse(T, numzeros), T, Ts, taumin, taumax, c, p, q);
    PilotsOnly<TruncatedSincPulse> ponly(P, pilots, TruncatedSincPulse(T, numzeros), T, Ts, taumin, taumax, c, p, q);
    dir.estimate(r);
    ponly.estimate(r);
    
    std::cout << "PilotOnly Ygrid test ... ";
    bool pass = true;
    vector<double> dygrid = dir.getYgrid();
    vector<double> pygrid = ponly.getYgrid();
    for(int i = 0; i < pygrid.size(); i++ ) {
        //std::cout << dygrid[i] << ", " << pygrid[i] << std::endl;
        pass &= fabs(dygrid[i] - pygrid[i]) < tol;
        pass &= fabs(pygrid[i] - std::abs(dir.Y(taumin + i*Delta))) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "PilotOnly SSgrid test ... ";
    pass = true;
    vector<double> pssgrid = ponly.getSSgrid();
    for(int i = 0; i < pygrid.size(); i++ ) {
        pass &= fabs(pssgrid[i] - ponly.SS(taumin + i*Delta)) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "PilotOnly Y test ... ";
    pass = true;
    for(double t = taumin; t < taumax; t+=0.1) {
        pass &= std::abs(ponly.Y(t) - dir.Y(t)) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
}

void test2() {
    Direct<TruncatedSincPulse> dir(P, D, pilots, TruncatedSincPulse(T, numzeros), T, Ts, taumin, taumax, c, p, q);
    TabulatedPilotsOnly<TruncatedSincPulse> ponly(P, pilots, TruncatedSincPulse(T, numzeros), T, Ts, taumin, taumax, c, p, q);
    dir.estimate(r);
    ponly.estimate(r);
    
    std::cout << "TabulatedPilotOnly Ygrid test ... ";
    bool pass = true;
    vector<double> dygrid = dir.getYgrid();
    vector<double> pygrid = ponly.getYgrid();
    for(int i = 0; i < pygrid.size(); i++ ) {
        //std::cout << dygrid[i] << ", " << pygrid[i] << std::endl;
        pass &= fabs(dygrid[i] - pygrid[i]) < tol;
        pass &= fabs(pygrid[i] - std::abs(dir.Y(taumin + i*Delta))) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "TabulatedPilotOnly SSgrid test ... ";
    pass = true;
    vector<double> pssgrid = ponly.getSSgrid();
    for(int i = 0; i < pygrid.size(); i++ ) {
        pass &= fabs(pssgrid[i] - ponly.SS(taumin + i*Delta)) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
    std::cout << "TabulatedPilotOnly Y test ... ";
    pass = true;
    for(double t = taumin; t < taumax; t+=0.1) {
        pass &= std::abs(ponly.Y(t) - dir.Y(t)) < tol;
    }
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
    
}


int main(int argc, char** argv) {
    
    for (int m = 0; m < numpilots; m++) P.push_back(m);
    for (int m = numpilots; m < L; m++) D.push_back(m);
    for (int m = 1; m <= L; m++) s.push_back(polar<double>(1.0, (2*pi*(rand()%M))/M));
    for ( int i = 0; i < P.size(); i++ ) pilots.push_back(s[P[i]]);
    for (int n = 1; n <= N; n++) r.push_back(polar<double>(1.0, -0.1 * n));
    
    test1();
    test2();

    return 1;
}

