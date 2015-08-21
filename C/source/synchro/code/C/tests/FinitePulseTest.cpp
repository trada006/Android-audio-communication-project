/* 
 * File:   FinitePulseTest.cpp
 * Author: Robby McKilliam
 *
 * Created on 23/01/2013, 2:19:26 PM
 */

#include <stdlib.h>
#include <iostream>
#include "FinitePulse.h"
#include "Util.h"

bool finiteSinc4zeros() {
    std::cout << "finiteSinc with 4 zeros ... ";
    double tol = 1e-9;
    double nzs = 4;
    double T = 1;
    TruncatedSincPulse sincpulse(T, nzs);
    bool pass = true;
    pass &= fabs(sincpulse.tmin() - -T * nzs) < tol;
    pass &= fabs(sincpulse.tmax() - T * nzs) < tol;
    pass &= fabs(sincpulse.pulse(-4.4) - complex<double>(0.0, 0.0)) < tol; //output from scala
    pass &= fabs(sincpulse.pulse(-3.9) - complex<double>(-0.025221324181627227, 0.0)) < tol;
    pass &= fabs(sincpulse.pulse(-3.4) - complex<double>(-0.08903843866360674, 0.0)) < tol;
    pass &= fabs(sincpulse.pulse(-2.9) - complex<double>(0.03391833252011936, 0.0)) < tol;
    pass &= fabs(sincpulse.pulse(-2.4) - complex<double>(0.12613778810677617, 0.0)) < tol;
    pass &= fabs(sincpulse.pulse(-1.9) - complex<double>(-0.05177008647807704, 0.0)) < tol;
    pass &= fabs(sincpulse.pulse(-1.9) - complex<double>(-0.05177008647807704, 0.0)) < tol;
    pass &= fabs(sincpulse.pulse(0.0) - complex<double>(1.0, 0.0)) < tol;
    pass &= fabs(sincpulse.pulse(4e-3) - complex<double>(0.9999973333354667, 0.0)) < tol;
    pass &= fabs(sincpulse.pulse(1.4) - complex<double>(-0.21623620818304484, 0.0)) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void normalisedTruncatedSinc() {
    std::cout << "normalised truncated sinc ... ";
    double tol = 1e-4;
    unsigned int nzs = 4;
    double T = 1;
    NormalisedFinitePulse<TruncatedSincPulse> p(TruncatedSincPulse(T, nzs));
    auto f = [&p](double t){return std::norm(p.pulse(t));};
    double energy1 = trapezoidal(f, p.tmin() - 0.3, p.tmax() + 0.3, 100000); //compute pulse energy
    double energy2 = trapezoidal(f, p.tmin(), p.tmax(), 100000); //compute pulse energy
    bool pass = true;
    pass &= fabs(1.0 - energy1) < tol;
    pass &= fabs(1.0 - energy2) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

bool finiteSinc2zeros() {
    std::cout << "finiteSinc with 2 zeros ... ";
    double tol = 1e-9;
    double nzs = 2;
    double T = 1;
    TruncatedSincPulse sincpulse(T, nzs);
    bool pass = true;
    pass &= fabs(sincpulse.tmin() - -T * nzs) < tol;
    pass &= fabs(sincpulse.tmax() - T * nzs) < tol;
    pass &= fabs(sincpulse.pulse(1.2) - complex<double>(-0.15591488063143982, 0.0)) < tol; //output from scala
    pass &= fabs(sincpulse.pulse(1.0) - complex<double>(0.0, 0.0)) < tol; //output from scala
    pass &= fabs(sincpulse.pulse(0.8) - complex<double>(0.23387232094715982, 0.0)) < tol; //output from scala
    pass &= fabs(sincpulse.pulse(0.6) - complex<double>(0.5045511524271047, 0.0)) < tol;
    pass &= fabs(sincpulse.pulse(0.4) - complex<double>(0.756826728640657, 0.0)) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void testRootRaisedCosineNearZero() {
    std::cout << "root raised cosine near zero... ";
    double tol = 1e-7;
    double nzs = 4;
    double T = 1.0;
    double beta = 0.5;
    TruncatedRootRaisedCosine p(T, beta, nzs);
    bool pass = true;
    pass &= fabs(1.136619772 - p.rootraisedcosine(0.0)) < tol; //from mathematica
    pass &= fabs(1.136617045 - p.rootraisedcosine(1e-3)) < tol; //from mathematica
    pass &= fabs(1.136617045 - p.rootraisedcosine(-1e-3)) < tol; //from mathematica
    pass &= fabs(1.136576129 - p.rootraisedcosine(-4e-3)) < tol; //from mathematica
    pass &= fabs(1.136576129 - p.rootraisedcosine(4e-3)) < tol; //from mathematica
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void testRootRaisedCosineNearBeta4() {
    std::cout << "root raised cosine near beta/4... ";
    double tol = 1e-6;
    double nzs = 4;
    double T = 1.0;
    double beta = 0.5;
    double t = 1.0 / 4 / beta;
    TruncatedRootRaisedCosine p(T, beta, nzs);
    bool pass = true;
    pass &= fabs(0.5786324696 - p.rootraisedcosine(t)) < tol; //from Mathematica
    pass &= fabs(0.5786324696 - p.rootraisedcosine(-t)) < tol; //from Mathematica

    pass &= fabs(0.5784538720 - p.rootraisedcosine(t + 1e-4)) < tol; //from Mathematica
    pass &= fabs(0.5788110635 - p.rootraisedcosine(t - 1e-4)) < tol; //from Mathematica
    pass &= fabs(0.5784538720 - p.rootraisedcosine(-t - 1e-4)) < tol; //from Mathematica
    pass &= fabs(0.5788110635 - p.rootraisedcosine(-t + 1e-4)) < tol; //from Mathematica

    pass &= fabs(0.5793468225 - p.rootraisedcosine(t - 4e-4)) < tol; //from Mathematica
    pass &= fabs(0.5779180564 - p.rootraisedcosine(t + 4e-4)) < tol; //from Mathematica
    pass &= fabs(0.5793468225 - p.rootraisedcosine(-t + 4e-4)) < tol; //from Mathematica
    pass &= fabs(0.5779180564 - p.rootraisedcosine(-t - 4e-4)) < tol; //from Mathematica
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void testRootRaisedCosine() {
    std::cout << "root raised cosine... ";
    double tol = 1e-6;
    double nzs = 4;
    double T = 1.0;
    double beta = 0.5;
    TruncatedRootRaisedCosine p(T, beta, nzs);
    std::vector<double> expected = {2.0 / 15 / pi, -1.0 / 3 / sqrt(2.0) / pi, -1.0 / 3 / pi, (pi + 2) / 2 / sqrt(2.0) / pi, 0.5 + 2 / pi,
        (pi + 2) / 2 / sqrt(2.0) / pi, -1.0 / 3 / pi, -1.0 / 3 / sqrt(2.0) / pi, 2.0 / 15 / pi}; //from Mathematica
    std::vector<double> test;
    for (double t = -2.0; t <= 2.0; t += 0.5) test.push_back(p.rootraisedcosine(t));
    bool pass = true;
    for (int i = 0; i < test.size(); i++)
        pass &= fabs(test[i] - expected[i]) < tol; //from Mathematica
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void testRootRaisedCosineFindZeros() {
    std::cout << "root raised cosine find zeros ... ";
    const double tol = 1e-4;
    const double T = 1.0;
    const double beta = 0.5;
    //test against roots found by Mathematica
    bool pass = true;
    pass &= fabs(0.873584 - TruncatedRootRaisedCosine(T, beta, 1).tmax()) < tol;
    pass &= fabs(1.69555 - TruncatedRootRaisedCosine(T, beta, 2).tmax()) < tol;
    pass &= fabs(-1.69555 - TruncatedRootRaisedCosine(T, beta, 2).tmin()) < tol;
    pass &= fabs(2.35734 - TruncatedRootRaisedCosine(T, beta, 3).tmax()) < tol;
    pass &= fabs(2.96409 - TruncatedRootRaisedCosine(T, beta, 4).tmax()) < tol;
    pass &= fabs(3.68054 - TruncatedRootRaisedCosine(T, beta, 5).tmax()) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

bool durationConstructedRRC() {
    std::cout << "construct rrc pulse with duration ... ";
    double tol = 1e-9;
    const double beta = 1.0/3.0;
    double duration = 10;
    int nzs = 7; //obtained by inspection.
    double T = 1;
    TruncatedRootRaisedCosine rrcpulse(T, beta, duration, true);
    TruncatedRootRaisedCosine rrcpulsetester(T, beta, nzs);
    bool pass = true;
    for(double t = -duration; t < duration; t+=0.01) 
        pass &= abs(rrcpulse.pulse(t) - rrcpulsetester.pulse(t)) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

int main(int argc, char** argv) {

    finiteSinc4zeros();
    finiteSinc2zeros();
    normalisedTruncatedSinc();

    testRootRaisedCosineNearZero();
    testRootRaisedCosineNearBeta4();
    testRootRaisedCosine();
    testRootRaisedCosineFindZeros();
    durationConstructedRRC();

    return (EXIT_SUCCESS);
}

