/* 
 * File:   OerderMeyerMasseyTest.cpp
 * Author: Robby McKilliam
 *
 * Created on 20/02/2013, 4:55:35 PM
 */

#include <stdlib.h>
#include <iostream>
#include <random>
#include <complex>
#include "FinitePulse.h"
#include "Util.h"
#include "Massey.h"
#include "FastAlgorithm.h"

/*
 * Simple C++ Test Suite
 */

using namespace std;


void testwithstore() {
  default_random_engine generator; 
  //generator.seed(10);
  const double benchtime = 10.0; //number of seconds to run each benchmark for.
  const double tol = 1e-2;

  const unsigned int M = 4;
  const unsigned int c = 4;
  const double T = 1.0; //symbol period
  const unsigned int p = 3;
  const unsigned int q = 17;
  const double Ts = p*T/q; //sample period
  const double taumin = 11.5;
  const double taumax = 110.0;
  
  const double tau0 = (taumax - taumin)/2 + taumin;
  const complex<double> a0 = polar<double>(1,0.1*pi);
  const int L = 100;
 

  //SNRs
  const double SNRdB = 20.0;
  const double SNR = pow(10.0, SNRdB/10.0);
  normal_distribution<double> noise(0.0,T/Ts/SNR/2);
  uniform_int_distribution<int> unifM(1, M); //for generating M-PSK symbols
  
   //const unsigned int numpilots = L/10; //about 10% pilots 
    const unsigned int numpilots = 30;
    //symbol position setup
    vector<int> P;
    for(int i = 0; i < numpilots; i++) P.push_back(i); //pilots at the front
    vector<int> D;
    for(int i = numpilots; i < L; i++) D.push_back(i); //data at the back
    vector<complex<double>> s;
    for (int m = 1; m <= L; m++) s.push_back(polar<double>(1.0, 2 * pi * unifM(generator) / M));
    vector<complex<double>> pilots;
    for (int m = 0; m < numpilots; m++) pilots.push_back(s[m]);

    TruncatedRootRaisedCosine tx(T,1.0/2.0,10);
    auto f = [] (complex<double> x) { return std::norm(x); }; //square law estimator
    Massey<decltype(tx)> est(P,D,pilots,tx,T,Ts,taumin,taumax,f);
    //Direct<decltype(tx)> est(P,D,pilots,tx,T,Ts,taumin,taumax,c,p,q);
  
        //number of samples
    unsigned int N = (unsigned int) ceil((T*(L+40)+taumax)/Ts); //number of samples
    
    //transmitted signal
    auto x = [&s,&tx,T,L] (double t) {
      int mini = max(1, (int)ceil((t - tx.tmax())/T));
      int maxi = min(L, (int)floor((t - tx.tmin())/T));
      complex<double> sum(0,0);
      for(int i = mini; i <= maxi; i++) sum += s[i-1] * tx.pulse(t - i*T);
      return sum;
    };
    
        //sampled received signal
    vector<complex<double>> r;
    for(int n = 1; n <= N; n++) r.push_back(a0*x(n*Ts-tau0) + complex<double>(noise(generator), noise(generator)));
    
    double tauhat = est.estimate(r);
    
    std::cout << "OerderMeyerAndMassey test ... ";
    bool pass = fabs(tau0 - tauhat) < tol;
    if (pass) std::cout << "PASS " << std::endl;
    else std::cout << "FAIL " << tau0 << ", " << tauhat << std::endl;
    
}

void testtabluated() {
  default_random_engine generator; 
  //generator.seed(10);
  const double benchtime = 10.0; //number of seconds to run each benchmark for.
  const double tol = 1e-2;

  const unsigned int M = 4;
  const unsigned int c = 4;
  const double T = 1.0; //symbol period
  const unsigned int p = 3;
  const unsigned int q = 17;
  const double Ts = p*T/q; //sample period
  const double taumin = 11.5;
  const double taumax = 110.0;
  
  const double tau0 = (taumax - taumin)/2 + taumin;
  const complex<double> a0 = polar<double>(1,0.1*pi);
  const int L = 100;
 

  //SNRs
  const double SNRdB = 20.0;
  const double SNR = pow(10.0, SNRdB/10.0);
  normal_distribution<double> noise(0.0,T/Ts/SNR/2);
  uniform_int_distribution<int> unifM(1, M); //for generating M-PSK symbols
  
   //const unsigned int numpilots = L/10; //about 10% pilots 
    const unsigned int numpilots = 30;
    //symbol position setup
    vector<int> P;
    for(int i = 0; i < numpilots; i++) P.push_back(i); //pilots at the front
    vector<int> D;
    for(int i = numpilots; i < L; i++) D.push_back(i); //data at the back
    vector<complex<double>> s;
    for (int m = 1; m <= L; m++) s.push_back(polar<double>(1.0, 2 * pi * unifM(generator) / M));
    vector<complex<double>> pilots;
    for (int m = 0; m < numpilots; m++) pilots.push_back(s[m]);

    TruncatedRootRaisedCosine tx(T,1.0/2.0,10);
    auto f = [] (complex<double> x) { return std::norm(x); }; //square law estimator
    TabulatedMassey<decltype(tx)> est(P,D,pilots,tx,T,Ts,taumin,taumax,f);
    //Direct<decltype(tx)> est(P,D,pilots,tx,T,Ts,taumin,taumax,c,p,q);
  
        //number of samples
    unsigned int N = (unsigned int) ceil((T*(L+40)+taumax)/Ts); //number of samples
    
    //transmitted signal
    auto x = [&s,&tx,T,L] (double t) {
      int mini = max(1, (int)ceil((t - tx.tmax())/T));
      int maxi = min(L, (int)floor((t - tx.tmin())/T));
      complex<double> sum(0,0);
      for(int i = mini; i <= maxi; i++) sum += s[i-1] * tx.pulse(t - i*T);
      return sum;
    };
    
        //sampled received signal
    vector<complex<double>> r;
    for(int n = 1; n <= N; n++) r.push_back(a0*x(n*Ts-tau0) + complex<double>(noise(generator), noise(generator)));
    
    double tauhat = est.estimate(r);
    
    std::cout << "Tabulated OerderMeyerAndMassey test ... ";
    bool pass = fabs(tau0 - tauhat) < tol;
    if (pass) std::cout << "PASS " << std::endl;
    else std::cout << "FAIL " << tau0 << ", " << tauhat << std::endl;
    
}

void testRecursivetabluated() {
  default_random_engine generator; 
  //generator.seed(10);
  const double benchtime = 10.0; //number of seconds to run each benchmark for.
  const double tol = 1e-2;

  const unsigned int M = 4;
  const unsigned int c = 4;
  const double T = 1.0; //symbol period
  const unsigned int p = 3;
  const unsigned int q = 17;
  const double Ts = p*T/q; //sample period
  const double taumin = 11.5;
  const double taumax = 110.0;
  
  const double tau0 = (taumax - taumin)/2 + taumin;
  const complex<double> a0 = polar<double>(1,0.1*pi);
  const int L = 100;
 

  //SNRs
  const double SNRdB = 20.0;
  const double SNR = pow(10.0, SNRdB/10.0);
  normal_distribution<double> noise(0.0,T/Ts/SNR/2);
  uniform_int_distribution<int> unifM(1, M); //for generating M-PSK symbols
  
   //const unsigned int numpilots = L/10; //about 10% pilots 
    const unsigned int numpilots = 30;
    //symbol position setup
    vector<int> P;
    for(int i = 0; i < numpilots; i++) P.push_back(i); //pilots at the front
    vector<int> D;
    for(int i = numpilots; i < L; i++) D.push_back(i); //data at the back
    vector<complex<double>> s;
    for (int m = 1; m <= L; m++) s.push_back(polar<double>(1.0, 2 * pi * unifM(generator) / M));
    vector<complex<double>> pilots;
    for (int m = 0; m < numpilots; m++) pilots.push_back(s[m]);

    TruncatedRootRaisedCosine tx(T,1.0/2.0,10);
    auto f = [] (complex<double> x) { return std::norm(x); }; //square law estimator
    RecursiveTabulatedMassey<decltype(tx)> est(P,D,pilots,tx,T,Ts,taumin,taumax,f);
    Massey<decltype(tx)> tester(P,D,pilots,tx,T,Ts,taumin,taumax,f);
  
        //number of samples
    unsigned int N = (unsigned int) ceil((T*(L+40)+taumax)/Ts); //number of samples
    
    //transmitted signal
    auto x = [&s,&tx,T,L] (double t) {
      int mini = max(1, (int)ceil((t - tx.tmax())/T));
      int maxi = min(L, (int)floor((t - tx.tmin())/T));
      complex<double> sum(0,0);
      for(int i = mini; i <= maxi; i++) sum += s[i-1] * tx.pulse(t - i*T);
      return sum;
    };
    
        //sampled received signal
    vector<complex<double>> r;
    for(int n = 1; n <= N; n++) r.push_back(a0*x(n*Ts-tau0) + complex<double>(noise(generator), noise(generator)));
    
    double tauhat = est.estimate(r);
    double tauhattester = tester.estimate(r);
    
    std::cout << "Recursive Tabulated OerderMeyerAndMassey test ... ";
    bool pass = fabs(tau0 - tauhat) < tol;
    pass = fabs(tauhattester - tauhat) < tol;
    if (pass) std::cout << "PASS " << std::endl;
    else std::cout << "FAIL " << tau0 << ", " << tauhat << std::endl;
    
}


int main(int argc, char** argv) {
    
    testwithstore();
    testtabluated();
    testRecursivetabluated();

    return (EXIT_SUCCESS);
}

