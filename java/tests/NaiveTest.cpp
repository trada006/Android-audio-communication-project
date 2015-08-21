/* 
 * File:   TimeOffsetEstimatorTest.cpp
 * Author: Robby McKilliam
 *
 * Created on 24/01/2013, 1:32:09 PM
 */

#include <cstdlib>
#include <iostream>

#include "TimeOffsetEstimator.h"
#include "FinitePulse.h"

static double tol = 1e-6;
static double T = 1.0; //symbol period
static unsigned int p = 1;
static unsigned int q = 5;
static unsigned int c = 5; //new c, its will divide T
static double Ts = p * T / q; //sample period
static double taumin = 10.0;
static double taumax = 30.0;
static double step = 5.0;
static unsigned int L = 10;
static unsigned int N = (int) ceil((T * L + taumax) / Ts);
static int numzeros = 2;
static vector<int> P;
static vector<int> D;
static vector<complex<double> > pilots;
static vector<complex<double> > r;
static unsigned int M = 4; //QPSK
static int numpilots = 4;

template <class Pulse> void NaiveDotr(Naive<Pulse>& est) {
    std::cout << "NaiveDotr ... ";
    bool pass = true;
    pass &= std::abs(est.m(taumin) - complex<double>(1.3499373963348864, 4.563483625803526)*Ts) < tol;
    pass &= std::abs(est.m(taumin+step) - complex<double>(1.6496241036561894, -4.463905699456121)*Ts) < tol;
    pass &= std::abs(est.m(taumin+2*step) - complex<double>(-3.9931090337278663, 2.5889754772421516)*Ts) < tol;
    pass &= std::abs(est.m(taumin+3*step) - complex<double>(4.748483513451538, 0.31562335065586683)*Ts) < tol;
    pass &= std::abs(est.m(taumin+4*step) - complex<double>(-3.6153254669352815, -3.0946947418331097)*Ts) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

template <class Pulse> void NaiveY(Naive<Pulse>& est) {
    std::cout << "NaiveY ... ";
    bool pass = true;
    pass &= std::abs(est.Y(taumin) - complex<double>(14.657852409140666, -3.228790691768866)*Ts) < tol;
    pass &= std::abs(est.Y(taumin+step) - complex<double>(-13.67538616337809, -6.185591310609769)*Ts) < tol;
    pass &= std::abs(est.Y(taumin+2*step) - complex<double>(7.254044220717808, 13.139884665524072)*Ts) < tol;
    pass &= std::abs(est.Y(taumin+3*step) - complex<double>(2.0523237347316847, -14.868278107005565)*Ts) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else {
        std::cout << "FAIL" << std::endl;
        std::cout << est.Y(taumin) << std::endl;
        std::cout << est.Y(taumin+step) << std::endl;
        std::cout << est.Y(taumin+2*step) << std::endl;
    }
}

template <class Pulse> void NaiveZ(Naive<Pulse>& est) {
    std::cout << "NaiveZ ... ";
    bool pass = true;
    pass &= fabs(est.Z(taumin) - 28.553768507361884*Ts) < tol;
    pass &= fabs(est.Z(taumin+step) - 28.553768507361887*Ts) < tol;
    pass &= fabs(est.Z(taumin+2*step) - 28.55376850736188*Ts) < tol;
    pass &= fabs(est.Z(taumin+3*step) - 28.553768507361884*Ts) < tol;
    pass &= fabs(est.Z(taumin+4*step) - 27.183844897660062*Ts) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else {
        std::cout << "FAIL" << std::endl;
        std::cout << est.Z(taumin) << std::endl;
        std::cout << est.Z(taumin+4*step) << std::endl;
    }
}

template <class Pulse> void NaiveCoarseEst(Naive<Pulse>& est) {
    std::cout << "Naive coarse estimate ... ";
    bool pass = fabs(29.199999999999932 - est.coarseMaximiseSS()) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else {
        std::cout << "FAIL" << std::endl;
        std::cout << est.coarseMaximiseSS() << std::endl;
    }
}

void NaiveFineEst(double tauhat) {
    std::cout << "Naive fine estimate ... ";
    bool pass = fabs(29.199999999999932 - tauhat) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else {
        std::cout << "FAIL" << std::endl;
        std::cout <<  tauhat << std::endl;
    }
}


int main(int argc, char** argv) {

    for (int m = 0; m < numpilots; m++) P.push_back(m);
    for (int m = numpilots; m < L; m++) D.push_back(m);
    for (int m = 1; m <= numpilots; m++) pilots.push_back(polar<double>(1.0, 0.1 * m));
    for (int n = 1; n <= N; n++) r.push_back(polar<double>(1.0, -0.1 * n));
    
    Naive<TruncatedSincPulse> est(P, D, pilots, TruncatedSincPulse(T, numzeros), T, Ts, taumin, taumax, c);
    double tauhat = est.estimate(r);
    
    NaiveDotr(est);
    NaiveY(est);
    NaiveZ(est);
    NaiveCoarseEst(est);
    NaiveFineEst(tauhat);
    
    return (EXIT_SUCCESS);
}

