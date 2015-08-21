/* 
 * File:   OerderMeyerMassey.h
 * Author: Robby McKilliam
 *
 * Created on 4 February 2013, 8:22 PM
 */

#ifndef OERDERMEYERMASSEY_H
#define	OERDERMEYERMASSEY_H

#include "TimeOffsetEstimator.h"
#include "TabulatedMainLoop.h"
#include "FinitePulse.h"
#include <complex>
#include <vector>
#include <algorithm>

/** 
 * Symbol timing estimator based on a nonlinear and filter, for example:
 * 
 * Oerder and Meyr. "Digital Filter and Square Timing Recovery" IEEE TRANSACTIONS ON COMMUNICATIONS, VOL 36. NO.5, MAY 1988
 * 
 * @params S is the indices of all the transmitted symbols 
 * @params g is the transmission pulse shape, outside the interval [gtmin, gtmax] the pulse g is zero.
 * @params Ts is the sample period
 * @params T is the symbol period
 * @params f_ is the nonlinear function to apply, should map complex to double.
 * @params N_ oversampling factor for matched filter.  Default is 4, which is suitable is most cases.
 */
template <class Pulse>
class NonlinearityAndFilter : public TimeOffsetEstimator {
public:

    NonlinearityAndFilter(const vector<int>& S,
            const Pulse& g,
            const double T,
            const double Ts,
            const double taumin,
            const double taumax,
            const function<double(complex<double>) > f_,
            const int N_ = 4) :
    g(g), T(T), Ts(Ts), taumin(taumin), taumax(taumax), f(f_), N(N_) {
        if(taumax <= taumin) throw "Maximum time offset taumax must be larger than minimum time offset taumin";
        Smax = *std::max_element(S.begin(), S.end());
        Smin = *std::min_element(S.begin(), S.end());
    }

    /**
     * Runs the Oerder Meyr estimator. This returns an estimator of the symbol period modulo T.
     */
    virtual double estimate(const vector<complex<double>>& r) {
        setupr(r);
        int kmin = (int) ceil(N * (g.tmin() + taumin + Smin * T) / T);
        int kmax = (int) floor(N * (g.tmax() + taumax + Smax * T) / T);
        //sum squared matched filtered output complex exponential
        complex<double> X(0, 0);
        for (int k = kmin; k <= kmax; k++)
	  X += polar(1.0, -(2 * pi * k) / N) * f(dotr(k * T / N));
        return -arg(X) / 2.0 / pi*T;
    }

    /** Inner product between the received signal and pulse g */
    virtual complex<double> dotr(double tau) const {
        int A = (int) ceil((g.tmin() + tau) / Ts);
        int B = (int) floor((g.tmax() + tau) / Ts);
        complex<double> sum = complex<double>(0, 0);
        for (int n = A; n <= B; n++) sum += r(n) * conj(g.pulse(n * Ts - tau));
        return sum;
    }

    virtual string name() const {
        return "NonlinearityAndFilter";
    }
    
    virtual void setupr(const vector<complex<double>>& r) {
        rmem = r.data();
        rmemsize = r.size();
    }

protected:
    int Smin;
    int Smax;
    const double T;
    const double Ts;
    const double taumin;
    const double taumax;
    const Pulse g;
    /** Pointer to current received signal */
    const complex<double>* rmem;
    /** size of current received signal */
    unsigned int rmemsize;
    /** Oversample factor on the initial matched filter.  Often this is 4. */
    const int N;
    /** The nonlinearity */
    const function<double(complex<double>)> f;
    
    /** The received sequence, zeros returned for unknown values */
    inline complex<double> r(int n) const {
        if (n > 0 && n <= rmemsize) return rmem[n - 1];
        else return complex<double>(0, 0);
    }

};

/**
 * Runs the Oerder Meyer estimator followed by Massey's frame synchronisation algorithm to allow for complex
 * symbols.  The original paper is:
 * 
 * JAMES L. MASSEY, "Optimum Frame Synchronization", IEEE TRANSACTIONS ON COMMUNICATIONS, VOL. COM-20, NO. 2, APRIL 1972.
 * 
 * This version does not store the rho_k sequence, a faster version is OerderMeyerAndMassey below
 */
template <class Pulse>
class MasseyNoStoreRho : public NonlinearityAndFilter<Pulse> {
    
public:

    MasseyNoStoreRho(const vector<int>& P,
            const vector<int>& D,
            const vector<complex<double> >& pilots,
            const Pulse& g,
            const double T,
            const double Ts,
            const double taumin,
            const double taumax,
            const function<double(complex<double>)> f,
            const int N = 4) : NonlinearityAndFilter<Pulse>(concatenate(P, D), g, T, Ts, taumin, taumax, f, N), P(P), D(D), pilots(pilots) {
    }

    virtual double estimate(const vector<complex<double>>& r) {
        setupgam(NonlinearityAndFilter<Pulse>::estimate(r)); //get estimate modulo the symbol period (call superclass)
        //int imin = (int) ceil((this->taumin - gam) / this->T);
        //int imax = (int) floor((this->taumax - gam) / this->T);
        int imin = (int) floor((this->taumin - gam) / this->T); //modified here, for some reason it seems to fail to get symbols at the boundary
        int imax = (int) ceil((this->taumax - gam) / this->T);
        int ibest = imin;
        double maxL = -1.0; //the objective is always positive, so this is fine.
        for (int i = imin; i <= imax; i++) {
            double thisL = L(i);
            if (thisL > maxL) {
                maxL = thisL;
                ibest = i;
            }
        }
        return ibest * this->T + gam;
    }

    /** The Massey-like objective function */
    double L(const int d) const {
        double dsum = 0.0;
        for (int i : D) dsum += abs(rho(i + d));
        complex<double> psum = complex<double>(0, 0);
        for (int i = 0; i < P.size(); i++) psum += rho(P[i] + d) * conj(pilots[i]);
        return dsum + abs(psum);
    }

protected:
    const vector<int> P;
    const vector<int> D;
    const vector<complex<double> > pilots;
    /** Variable contains the symbol time estimate.  MUTABLE */
    double gam;
    
    virtual void setupgam(double gam) {
        this->gam = gam;
    }
    
    virtual complex<double> rho(const int k) const {
        return this->dotr((k + 1) * this->T + gam);
    }

};


/**
 * Runs the Oerder Meyer estimator followed by Massey's frame synchronisation algorithm to allow for complex
 * symbols.  The original paper is:
 * 
 * JAMES L. MASSEY, "Optimum Frame Synchronization", IEEE TRANSACTIONS ON COMMUNICATIONS, VOL. COM-20, NO. 2, APRIL 1972.
 * 
 * This version does not store the rho_k sequence, a faster version is OerderMeyerAndMassey below
 */
template <class Pulse>
class Massey : public MasseyNoStoreRho<Pulse> {
    
public:

    Massey(const vector<int>& P,
            const vector<int>& D,
            const vector<complex<double>>& pilots,
            const Pulse& g,
            const double T,
            const double Ts,
            const double taumin,
            const double taumax,
            const function<double(complex<double>) > f,
            const int N = 4) : MasseyNoStoreRho<Pulse>(P,D,pilots,g,T,Ts,taumin,taumax,f,N) {
        int imin = (int) floor((this->taumin - 0.5) / this->T); //modified here, for some reason it seems to fail to get symbols at the boundary
        int imax = (int) ceil((this->taumax + 0.5) / this->T);
        int kmin = imin + this->Smin;
        int kmax = imax + this->Smax;
        rhostore.resize(kmax - kmin + 1); //rhostore will never need to be bigger that this
    }

    virtual string name() const {
        return "Massey";
    }

protected:
    vector<complex<double>> rhostore;
    int rhooffset;
    
    /** Override setgam to also store all the values of rho we will need */
    virtual void setupgam(double gam) {
        this->gam = gam;
        //int imin = (int) ceil((this->taumin - gam) / this->T);
        //int imax = (int) floor((this->taumax - gam) / this->T);
        int imin = (int) floor((this->taumin - this->gam) / this->T); //modified here, for some reason it seems to fail to get symbols at the boundary
        int imax = (int) ceil((this->taumax - this->gam) / this->T);
        int kmin = imin + this->Smin;
        int kmax = imax + this->Smax;
        rhooffset = kmin;
        for(int k = kmin; k <= kmax; k++) 
            rhostore[k - rhooffset] = MasseyNoStoreRho<Pulse>::rho(k); //fill rhostore
    }
    
    /** Override rho to take from the array of stored values */
    inline virtual complex<double> rho(const int k) const {
        return rhostore[k - rhooffset];
    }
            
};

#define DEFAULTTABLESIZE 100000

/**
 * Tabulated Oerder Meyer.  Store pulse in a lookup table for fast access.  Table is memory
 * aligned for extra speed.
 */
template <class Pulse>
class TabulatedMassey : public Massey<Pulse> {
    
public:

   TabulatedMassey(const vector<int>& P,
            const vector<int>& D,
            const vector<complex<double>>& pilots,
            const Pulse& g,
            const double T,
            const double Ts,
            const double taumin,
            const double taumax,
            const function<double(complex<double>) > f,
            const int N = 4,
            const unsigned int mintabs = DEFAULTTABLESIZE) : 
    Massey<Pulse>(P,D,pilots,g,T,Ts,taumin,taumax,f,N) {
        tabls = (int) ceil(mintabs*Ts/T);
        stepwidth = Ts/tabls;
        //fill up the pulse table with the transmit pulse
        for (double t = g.tmin(); t <= g.tmax(); t += stepwidth) pulsetable.push_back(conj(g.pulse(t)));
    }
    
    /** Inner product between the received signal and pulse g */
    virtual complex<double> dotr(double tau) const {
        //setup loop counters
        int A = (int) ceil((this->g.tmin() + tau) / this->Ts);
        int B = (int) floor((this->g.tmax() + tau) / this->Ts);
        int nfrom = max(1, A);
        int nto = min((int)this->rmemsize,B);
        double startt = nfrom * this->Ts - tau;
        int ifrom = (int)round(( startt - this->g.tmin() ) / stepwidth);
        int istep = tabls;
        int mfrom = nfrom - 1;
        int mto = nto - 1;
        int mstep = 1;
        
        return mainFilterLoop(mfrom,mto,mstep,ifrom,istep,this->rmem,&pulsetable[0]);
        
    }

    virtual string name() const {
        return "TabulatedMassey";
    }

protected:
    /** Table stores values of the transmit pulse g.pulse */
    vector<complex<double>> pulsetable;
    /**  Width in t between elements in the table. */
    double stepwidth;
    /** Multiplier for indexing the table */
    int tabls;
};

/**
 * Recursive Tabulated Oerder Meyer.  Exploits the fact that data symbols are in a contiguous
 * block similarly to Recursive
 */
template <class Pulse>
class RecursiveTabulatedMassey : public TabulatedMassey<Pulse> {
    
public:

    RecursiveTabulatedMassey(const vector<int>& P,
            const vector<int>& D,
            const vector<complex<double> >& pilots,
            const Pulse& g,
            const double T,
            const double Ts,
            const double taumin,
            const double taumax,
            const function<double(complex<double>) > f,
            const int N = 4,
            const unsigned int mintabs = DEFAULTTABLESIZE) : 
    TabulatedMassey<Pulse>(P,D,pilots,g,T,Ts,taumin,taumax,f,N,mintabs),
        Dmin(*std::min_element(D.begin(), D.end())),
        Dmax(*std::max_element(D.begin(), D.end())) {
        //check D is contiguous and sorted
        for (unsigned int i = 0; i < D.size() - 1; i++) {
            if (D[i] != (D[i + 1] - 1)) {
                std::ostringstream str;
                str << "Exception constructing RecursiveTabulatedOerderMeyerAndMassey. The data symbols must be contiguous." << endl;
                throw str.str();
            }
        }
        int imin = (int) floor((taumin - 0.5) / T);
        int imax = (int) ceil((taumax + 0.5) / T);
        dgrid.resize(imax - imin + 1, 0);
    }
    
    virtual double estimate(const vector<complex<double> >& r) {
        this->setupgam(NonlinearityAndFilter<Pulse>::estimate(r)); //get estimate modulo the symbol period (call superclass)
        filldgrid();
        int imin = (int) floor((this->taumin - this->gam) / this->T);
        int imax = (int) ceil((this->taumax - this->gam) / this->T);
        int ibest = imin;
        double maxL = -1; //the objective is always positive, so this is fine.
        for (int i = imin; i <= imax; i++) {
            double thisL = this->L(i);
            if (thisL > maxL) {
                maxL = thisL;
                ibest = i;
            }
        }
        return ibest * this->T + this->gam;
    }

    /** The Massey-like objective function */
    virtual double L(const int d) const {
        std::complex<double> psum = std::complex<double>(0, 0);
        for (int i = 0; i < this->P.size(); i++) psum += this->rho(this->P[i] + d) * std::conj<double>(this->pilots[i]);
        return dfunc(d) + std::abs<double>(psum);
    }
    
    void filldgrid() {
        imin = (int) floor((this->taumin - this->gam) / this->T);
        int imax = (int) ceil((this->taumax - this->gam) / this->T);
        int dto = imax - imin;
        for (int i  = 0; i < this->D.size(); i++) dgrid[0] += std::abs<double>(this->rho(this->D[i] + imin));
        for(int k = 0; k < dto; k++)
            dgrid[k+1] = dgrid[k] + std::abs<double>(this->rho(Dmax+imin+k+1)) - std::abs<double>(this->rho(Dmin+imin+k));    
    }
    
    inline double dfunc(int i) const { 
        return dgrid[i - imin]; 
    }

    virtual string name() const {
        return "RecursiveTabulatedOerderMeyrAndMassey";
    }
    
    protected:
        vector<double> dgrid;
        int imin;
        /** Minimum data symbol index */
        const int Dmin;
        /** Maximum data symbol index */
        const int Dmax;
            
};

#endif	/* OERDERMEYERMASSEY_H */

