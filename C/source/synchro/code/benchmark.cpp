/* 
 * File:   benchmark.cpp
 * Author: Robby McKilliam
 */

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "Util.h"
#include "FinitePulse.h"
#include "TimeOffsetEstimator.h"
#include "FastAlgorithm.h"
#include "TabulatedFastAlgorithm.h"
#include "Massey.h"
#include <random>
#include <ctime>
#include <functional>

using namespace std;


  static constexpr double benchtime = 10.0; //number of seconds to run each benchmark for.
  static constexpr unsigned int M = 4;
  static constexpr unsigned int c = 15;
  static constexpr double T = 1.0; //symbol period
  static constexpr unsigned int p = 2;
  static constexpr unsigned int q = 11;
  static constexpr double Ts = p*T/q; //sample period
  static constexpr double taumin = 50*T;
  static constexpr double taumax = 150*T;
static constexpr double tau0 = 82.1111*T;

/* const unsigned int M = 4;
  const unsigned int c = 12;
  const VALTYPE T = 62.5e-6; //symbol period
  const unsigned int p = 1;
  const unsigned int q = 4;
  const VALTYPE Ts = p*T/q; //sample period
  //const VALTYPE taumin = 10.0;
  //const VALTYPE taumax = 110.0;
  //Calculation of taumin and taumax (as defined in Matlab)
  const VALTYPE OS = 4 ; //Oversampling Factor
  const VALTYPE Symbol_Rate = 64000; 
  const VALTYPE LengthRRCCoefficient = 29;
  const VALTYPE TauMin = 0;
  const VALTYPE TauMax = 166.5*T;
  const VALTYPE PulseDuration = (LengthRRCCoefficient - 1)/Symbol_Rate;
  const VALTYPE TauCorrection = (PulseDuration/2) - (OS-1)*T/OS;
  //TauMin and TauMax as defined in Matlab
  const VALTYPE taumin = TauMin + TauCorrection;
  const VALTYPE taumax = TauMax + TauCorrection;
  const VALTYPE tau0 = (taumax - taumin)/2 + taumin;*/

template <class Pulse, class Estimator> 
void runbenchmark(const vector<int>& Ls, const Pulse& tx, function<Estimator(vector<int>&,vector<int>&,vector<complex<double>>&)> estf, string name) {

  default_random_engine generator; 
      //SNRs
  double SNRdB = 10.0;
  double SNR = pow(10.0, SNRdB/10.0);
  normal_distribution<double> noise(0.0,T/Ts/SNR/2);
  uniform_int_distribution<int> unifM(1, M); //for generating M-PSK symbols
  
  vector<double> times; //array to store output
  
  for(auto L : Ls) {

    const unsigned int numpilots = L/5; //about 10% pilots 
    //symbol position setup
    vector<int> P;
    for(int i = 0; i < numpilots; i++) P.push_back(i); //pilots at the front
    vector<int> D;
    for(int i = numpilots; i < L; i++) D.push_back(i); //data at the back
    vector<complex<double>> s;
    for (int m = 1; m <= L; m++) s.push_back(polar<double>(1.0, 2 * pi * unifM(generator) / M));
    vector<complex<double>> pilots;
    for (int m = 0; m < numpilots; m++) pilots.push_back(s[m]);

    //construct the estimator we will run
    Estimator est = estf(P,D,pilots);
   
    //number of samples
    unsigned int N = (unsigned int) ceil((T*(L+40)+taumax)/Ts); //number of samples
    const complex<double> a0 = polar<double>(1,0.1*pi); //phase and amplitude
    
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

    cout << "Benchmarking " << name << " L = " << L << " ... ";
    long iters = 0;
    double errmse = 0.0;
    clock_t started = clock();
    while( ((double)(clock() - started))/CLOCKS_PER_SEC < benchtime) {
      double tauhat = est.estimate(r);
      errmse += (tauhat - tau0)*(tauhat - tau0); //use tauhat otherwise it might get compiled out!
      iters++;
    }
    clock_t stopped = clock();
    double microsecs = ((double)(stopped - started))/CLOCKS_PER_SEC/iters/L*1000000;
    times.push_back(microsecs);
    cout << " requires  " << microsecs << " microseconds per symbol with average error " << (errmse/iters) << endl;
  }

  ofstream file(string("data/") + name + string("bench"));
  for(int i = 0; i < Ls.size(); i++) file << Ls[i] << "\t" << times[i] << endl;
  file.close();

}

int main(int argc, char** argv) {

  //const vector<int> Ls = {20,30,40,50,60,80,100,140,200,250,300,400,500,600,800,1000,1200,1500,2000,2600,3250,4000};
  const vector<int> Ls = {50,60,80,100,140,200,250,300,400,500,600,800,1000,1200,1500,2000,2600,3250,4000};
  //const vector<int> Ls = {11,50,1000,4000};
  //const vector<int> Ls = {3990};

  //the transmit pulse
  TruncatedRootRaisedCosine tx(T,1.0/3.0,15);

  //functions that construct estimators
  auto rec = [&] (vector<int>& P, vector<int>& D, vector<complex<double>>& pilots) 
    { return TabulatedRecursive<TruncatedRootRaisedCosine>(P,D,pilots,tx,T,Ts,taumin,taumax,c,p,q); };
  auto dir = [&] (vector<int>& P, vector<int>& D, vector<complex<double>>& pilots) 
    { return TabulatedDirect<decltype(tx)>(P,D,pilots,tx,T,Ts,taumin,taumax,c,p,q); };
  auto recom = [&] (vector<int>& P, vector<int>& D, vector<complex<double>>& pilots) 
    { return RecursiveTabulatedMassey<decltype(tx)>(P,D,pilots,tx,T,Ts,taumin,taumax, [] (complex<double> x) { return std::norm(x); }); };
  auto recnorefine = [&] (vector<int>& P, vector<int>& D, vector<complex<double>>& pilots) 
    { return TabulatedRecursiveNoRefine<TruncatedRootRaisedCosine>(P,D,pilots,tx,T,Ts,taumin,taumax,c,p,q); };

  //run the benchmarks
  runbenchmark<decltype(tx),TabulatedRecursive<decltype(tx)>>(Ls,tx,rec,"Recursive");
  runbenchmark<decltype(tx),RecursiveTabulatedMassey<decltype(tx)>>(Ls,tx,recom,"RecursiveTabulatedMassey");
  runbenchmark<decltype(tx),TabulatedRecursiveNoRefine<decltype(tx)>>(Ls,tx,recnorefine,"RecursiveNoRefinement");
  
}
