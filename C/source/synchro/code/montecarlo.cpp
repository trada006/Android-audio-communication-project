/* 
 * File:   montecarlo.cpp
 * Author: Robby McKilliam
 */

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "Util.h"
#include "FinitePulse.h"
#include "FastAlgorithm.h"
#include "TabulatedFastAlgorithm.h"
#include "Massey.h"
#include "TimeOffsetEstimator.h"
#include <random>
#include <functional>
#include <future>
#include <chrono>

using namespace std;

typedef std::chrono::duration<double, std::ratio<1>> seconds;

static constexpr unsigned int M = 4; //QPSK
static constexpr unsigned int iters = 50;
static constexpr double T = 1.0;
static constexpr unsigned int p = 2;
static constexpr unsigned int q = 11;
static constexpr double Ts = p*T/q;
static constexpr double taumin = 50.0;
static constexpr double taumax = 150.0;
static constexpr unsigned int c = 15;

template <class Pulse, class Estimator> 
void runsim(const unsigned int L, const int mseqN, const vector<double>& snrdbs, const Pulse& tx, function<Estimator(vector<int>&,vector<int>&,vector<complex<double>>&,const Pulse&)> estf, const string name) {

  cout << "Running simulation " << name << " with L = " << L << " " << flush;
  auto starttime = chrono::system_clock::now();

  //a list of noise distributions we will use
  vector<normal_distribution<double>> noises;
  for(double snrdb : snrdbs) {
    double snr = pow(10.0, snrdb/10.0);
    noises.push_back(normal_distribution<double>(0.0,sqrt(T/Ts/snr/2)));
  }

  vector<future<double>> msearray(snrdbs.size()); //store the mses
  //vector<double> msearray(snrdbs.size()); //store the mses
  for( int snrind = 0; snrind < snrdbs.size(); snrind++ ) {
    //run each snr in a separate thread
    msearray[snrind] = async( launch::async, [=] {
	auto noise = noises[snrind];
	default_random_engine generator((long)clock()); //random number generator, seed with clock 
	uniform_int_distribution<int> unifM(1, M); //for generating M-PSK symbols
	uniform_real_distribution<double> unif(taumin,taumax); //generator for true timeoffset 
	
	const unsigned int N = (unsigned int) ceil((T*(L+40)+taumax)/Ts); //number of samples
	const complex<double> a0 = polar<double>(1,0.1*pi); //phase and amplitude

	vector<complex<double>> r(N); //sampled received signal
	vector<complex<double>> s(L); //stores transmitted signal
	vector<complex<double>> pilots = pilotmsequence(mseqN); //stores array of pilots
	int numpilots = pilots.size();
	for(int i = 0; i < numpilots; i++) s[i] = pilots[i]; //first symbols with pilots
    
	//symbol position setup
	vector<int> P;
	for(int i = 0; i < numpilots; i++) P.push_back(i); //pilots at the front
	vector<int> D;
	for(int i = numpilots; i < L; i++) D.push_back(i); //data at the back

	//function generates transmitted signal
	auto x = [&s,&tx,L] (double t) {
	  int mini = max((int)1, (int)ceil((t - tx.tmax())/T));
	  int maxi = min((int)L, (int)floor((t - tx.tmin())/T));
	  complex<double> sum(0,0);
	  for(int i = mini; i <= maxi; i++) sum += s[i-1] * tx.pulse(t - i*T); //this here is already fairly time consuming!
	  return sum;
	};
    
	Estimator est = estf(P,D,pilots,tx); //construct the estimator we will run
    
	double mse = 0.0;
	for(int itr = 0; itr < iters; itr++) {
	  double tau0 = unif(generator);
	  for(int i = numpilots; i < L; i++) s[i] = polar<double>(1.0, 2 * pi * unifM(generator) / M); //fill s, random M-PSK
	  for(int n = 1; n <= N; n++) r[n-1] = a0*x(n*Ts-tau0) + complex<double>(noise(generator), noise(generator)); //generate received signal
	  double tauhat = est.estimate(r);
	  mse+=(tauhat - tau0)*(tauhat - tau0);
	}
	cout << "." << flush;	
	return mse/iters;
	} );
    //msearray[snrind] = mse/iters;
  }

  //now write to file
  for(int i = 0; i < msearray.size(); i++) msearray[i].wait(); //wait for all the threads to finish, not strictly necessary but avoid overwriting files.
  ofstream file(string("data/") + name);      
  for( int snrind = 0; snrind < snrdbs.size(); snrind++ ) 
    file << snrdbs[snrind] << "\t" << msearray[snrind].get() << endl;
  file.close();

  auto endtime = chrono::system_clock::now();
  auto elapsed_seconds = std::chrono::duration_cast<seconds>(endtime-starttime).count();
  cout << " finished in  " << elapsed_seconds << " seconds." << endl;

}


template <class Pulse> 
void runAllPilotsSim(const unsigned int L, const int mseqN, const vector<double>& snrdbs, const Pulse& tx, const string name) {

  cout << "Running simulation " << name << " with L = " << L << " " << flush;
  auto starttime = chrono::system_clock::now();

  //a list of noise distributions we will use
  vector<normal_distribution<double>> noises;
  for(double snrdb : snrdbs) {
    double snr = pow(10.0, snrdb/10.0);
    noises.push_back(normal_distribution<double>(0.0,sqrt(T/Ts/snr/2)));
  }

  vector<future<double>> msearray(snrdbs.size()); //store the mses
  //vector<double> msearray(snrdbs.size()); //store the mses
  for( int snrind = 0; snrind < snrdbs.size(); snrind++ ) {
    //run each snr in a separate thread
    msearray[snrind] = async( launch::async, [=] {
	auto noise = noises[snrind];
	default_random_engine generator((long)clock()); //random number generator, seed with clock 
	uniform_int_distribution<int> unifM(1, M); //for generating M-PSK symbols
	uniform_real_distribution<double> unif(taumin,taumax); //generator for true timeoffset 
	
	const unsigned int N = (unsigned int) ceil((T*(L+40)+taumax)/Ts); //number of samples
	const complex<double> a0 = polar<double>(1,0.1*pi); //phase and amplitude

	vector<complex<double>> r(N); //sampled received signal
	vector<complex<double>> s(L); //stores transmitted signal
	vector<complex<double>> pilots = pilotmsequence(mseqN); //stores array of pilots 
	int numpilots = pilots.size();
	for(int i = 0; i < numpilots; i++) s[i] = pilots[i]; //first symbols with pilots
    
	//symbol position setup
	vector<int> P;
	for(int i = 0; i < L; i++) P.push_back(i); //every symbol is a pilot

	//function generates transmitted signal
	auto x = [&s,&tx,L] (double t) {
	  int mini = max((int)1, (int)ceil((t - tx.tmax())/T));
	  int maxi = min((int)L, (int)floor((t - tx.tmin())/T));
	  complex<double> sum(0,0);
	  for(int i = mini; i <= maxi; i++) sum += s[i-1] * tx.pulse(t - i*T); //this here is already fairly time consuming!
	  return sum;
	};

	TabulatedPilotsOnly<Pulse> est(P,s,tx,T,Ts,taumin,taumax,c,p,q); //construct estimator where every symbol is a pilot
    
	double mse = 0.0;
	for(int itr = 0; itr < iters; itr++) {
	  for(int i = numpilots; i < L; i++) s[i] = polar<double>(1.0, 2 * pi * unifM(generator) / M); //fill s, random M-PSK
	  double tau0 = unif(generator);
	  for(int n = 1; n <= N; n++) r[n-1] = a0*x(n*Ts-tau0) + complex<double>(noise(generator), noise(generator)); //generate received signal
	  est.setPilots(s); //set the data symbols as pilots
	  double tauhat = est.estimate(r);
	  mse+=(tauhat - tau0)*(tauhat - tau0);
	}
	cout << "." << flush;	
	return mse/iters;
      } );
    //msearray[snrind] = mse/iters;
  }

  //now write to file
  for(int i = 0; i < msearray.size(); i++) msearray[i].wait(); //wait for all the threads to finish, not strictly necessary but avoid overwriting files.
  ofstream file(string("data/") + name);      
  for( int snrind = 0; snrind < snrdbs.size(); snrind++ ) 
    file << snrdbs[snrind] << "\t" << msearray[snrind].get() << endl;
  file.close();

  auto endtime = chrono::system_clock::now();
  auto elapsed_seconds = std::chrono::duration_cast<seconds>(endtime-starttime).count();
  cout << " finished in  " << elapsed_seconds << " seconds." << endl;

}

typedef NormalisedFinitePulse<TruncatedRootRaisedCosine> nrrc;

int main(int argc, char** argv) {
    
    vector<double> snrdbs;//snrs in db we will run
    for(int db = -26; db <= 22; db+=2) snrdbs.push_back(db);
  
    //estimators we will benchmark
    auto sln = [&] (vector<int>& P, vector<int>& D, vector<complex<double>>& pilots, const nrrc& tx) //square then Massey
      { return TabulatedMassey<nrrc>(P,D,pilots,tx,T,Ts,taumin,taumax,[] (complex<double> x) { return std::norm(x); }); };
    auto avn = [&] (vector<int>& P, vector<int>& D, vector<complex<double>>& pilots, const nrrc& tx) //absolute value then Massey
      { return TabulatedMassey<nrrc>(P,D,pilots,tx,T,Ts,taumin,taumax,[] (complex<double> x) { return std::abs(x); }); };
    auto rec = [&] (vector<int>& P, vector<int>& D, vector<complex<double>>& pilots, const nrrc& tx) //pilot and data (new)
      { return TabulatedRecursive<nrrc>(P,D,pilots,tx,T,Ts,taumin,taumax,c,p,q); };
    auto pilotonly = [&] (vector<int>& P, vector<int>& D, vector<complex<double>>& pilots, const nrrc& tx) //pilot only correlator
      { return TabulatedPilotsOnly<nrrc>(P,pilots,tx,T,Ts,taumin,taumax,c,p,q); };
     
    { //run large rolloff
      nrrc tx(TruncatedRootRaisedCosine(T,1.0/3.0,15));  //the transmit pulse
      runsim<decltype(tx), TabulatedRecursive<decltype(tx)>>(75, 4, snrdbs, tx, rec, "recursiveP15L75");
      runsim<decltype(tx), TabulatedRecursive<decltype(tx)>>(2555, 9, snrdbs, tx, rec, "recursiveP511L2555");
      //runsim<decltype(tx), TabulatedRecursive<decltype(tx)>>(5115, 10, snrdbs, tx, rec, "recursiveP1023L5115");
      //runsim<decltype(tx), TabulatedRecursive<decltype(tx)>>(10235, 11, snrdbs, tx, rec, "recursiveP2047L10235");
      runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(75, 4, snrdbs, tx, sln, "slnP15L75");
      runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(2555, 9, snrdbs, tx, sln, "slnP511L2555");
      //runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(5115, 10, snrdbs, tx, sln, "slnP1023L5115");
      //runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(10235, 11, snrdbs, tx, sln, "slnP2047L10235");
      runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(75, 4, snrdbs, tx, avn, "avnP15L75");
      runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(2555, 9, snrdbs, tx, avn, "avnP511L2555");
      //runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(5115, 10, snrdbs, tx, avn, "avnP1023L5115");
      //runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(10235, 11, snrdbs, tx, avn, "avnP2047L10235");
      runsim<decltype(tx), TabulatedPilotsOnly<decltype(tx)>>(75, 4, snrdbs, tx, pilotonly, "pilotonlyP15L75");
      runsim<decltype(tx), TabulatedPilotsOnly<decltype(tx)>>(2555, 9, snrdbs, tx, pilotonly, "pilotonlyP511L2555");
      //runsim<decltype(tx), TabulatedPilotsOnly<decltype(tx)>>(5115, 10, snrdbs, tx, pilotonly, "pilotonlyP1023L5115");
      //runsim<decltype(tx), TabulatedPilotsOnly<decltype(tx)>>(10235, 11, snrdbs, tx, pilotonly, "pilotonlyP2047L10235");
      runAllPilotsSim<decltype(tx)>(75, 4, snrdbs, tx, "allpilotsP15L75");
      runAllPilotsSim<decltype(tx)>(2555, 9, snrdbs, tx, "allpilotsP511L2555");
      //runAllPilotsSim<decltype(tx)>(5115, 10, snrdbs, tx, "allpilotsP1023L5115");
      //runAllPilotsSim<decltype(tx)>(10235, 11, snrdbs, tx, "allpilotsP2047L10235");
    }
    
    { //run small rolloff
      nrrc tx(TruncatedRootRaisedCosine(T,1.0/30,30));  //the transmit pulse
      runsim<decltype(tx), TabulatedRecursive<decltype(tx)>>(75, 4, snrdbs, tx, rec, "recursiveP15L75smallrolloff");
      runsim<decltype(tx), TabulatedRecursive<decltype(tx)>>(2555, 9, snrdbs, tx, rec, "recursiveP511L2555smallrolloff");
      //runsim<decltype(tx), TabulatedRecursive<decltype(tx)>>(5115, 10, snrdbs, tx, rec, "recursiveP1023L5115smallrolloff");
      //runsim<decltype(tx), TabulatedRecursive<decltype(tx)>>(10235, 11, snrdbs, tx, rec, "recursiveP2047L10235smallrolloff");
      runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(75, 4, snrdbs, tx, sln, "slnP15L75smallrolloff");
      runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(2555, 9, snrdbs, tx, sln, "slnP511L2555smallrolloff");
      //runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(5115, 10, snrdbs, tx, sln, "slnP1023L5115smallrolloff");
      //runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(10235, 11, snrdbs, tx, sln, "slnP2047L10235smallrolloff");
      runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(75, 4, snrdbs, tx, avn, "avnP15L75smallrolloff");
      runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(2555, 9, snrdbs, tx, avn, "avnP511L2555smallrolloff");
      //runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(5115, 10, snrdbs, tx, avn, "avnP1023L5115smallrolloff");
      //runsim<decltype(tx), TabulatedMassey<decltype(tx)>>(10235, 11, snrdbs, tx, avn, "avnP2047L10235smallrolloff");
      runsim<decltype(tx), TabulatedPilotsOnly<decltype(tx)>>(75, 4, snrdbs, tx, pilotonly, "pilotonlyP15L75smallrolloff");
      runsim<decltype(tx), TabulatedPilotsOnly<decltype(tx)>>(2555, 9, snrdbs, tx, pilotonly, "pilotonlyP511L2555smallrolloff");
      //runsim<decltype(tx), TabulatedPilotsOnly<decltype(tx)>>(5115, 10, snrdbs, tx, pilotonly, "pilotonlyP1023L5115smallrolloff");
      //runsim<decltype(tx), TabulatedPilotsOnly<decltype(tx)>>(10235, 11, snrdbs, tx, pilotonly, "pilotonlyP2047L10235smallrolloff");
      runAllPilotsSim<decltype(tx)>(75, 4, snrdbs, tx, "allpilotsP15L75smallrolloff");
      runAllPilotsSim<decltype(tx)>(2555, 9, snrdbs, tx, "allpilotsP511L2555smallrolloff");
      //runAllPilotsSim<decltype(tx)>(5115, 10, snrdbs, tx, "allpilotsP1023L5115smallrolloff");
      //runAllPilotsSim<decltype(tx)>(10235, 11, snrdbs, tx, "allpilotsP2047L10235smallrolloff");
    }
    
}
