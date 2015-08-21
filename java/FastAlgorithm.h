/* 
 * File:   FastAlgorithm.h
 * Author: Robby McKilliam
 *
 * Created on 4 February 2013, 8:13 PM
 */

#ifndef FASTALGORITHM_H
#define	FASTALGORITHM_H

#include "TimeOffsetEstimator.h"
#include "FinitePulse.h"
#include <complex>
#include <vector>
#include <algorithm>

/**
 * Override Naive with a faster way to compute the v vectors.  This does not store the v_k sequence.
 */
template <class Pulse> //P should be a FinitePulse
class FastAlgorithmNoStorev : public Naive<Pulse> {
public:

    FastAlgorithmNoStorev(const vector<int>& P,
            const vector<int>& D,
            const vector<complex<double> >& pilots,
            const Pulse& g,
            const double T,
            const double Ts,
            const double taumin,
            const double taumax,
            const unsigned int c,
            const unsigned int p,
            const unsigned int q,
            const double brenttol = DEFAULTBRENTTOL) :
    Naive<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, brenttol),
    p(p),
    q(q),
    K(floor((taumax - taumin) / this->Delta)),
    Zgrid(K, 0.0),
    Ygridmag(K, 0.0),
    SSgrid(K, 0.0) {
        //setup a, b, n0 and m0 for the polyphase filter
        int d = gcd(q, c * p);
        a = c * p / d;
        b = q / d;
        extended_gcd(b, a, d, n0, m0);
    }

    virtual double coarseMaximiseSS() {
        fillSSgrid(); //compute objective function on the grid
        auto biggest = std::max_element(SSgrid.begin(), SSgrid.end()); //get iterator to max element
        int khat = std::distance(SSgrid.begin(), biggest); //get position of max element
        return this->taumin + khat*this->Delta;
    }

    /** fills the vector SSgrid with values of the objective function SS */
    void fillSSgrid() {
        fillZgrid();
        fillYgrid();
        for (unsigned int i = 0; i < K; i++) SSgrid[i] = Zgrid[i] + Ygridmag[i];
    }

    /** zero filled received signal z_n in the paper */
    inline complex<double> zn(int n) const {
        if (n % a == 0) return this->r(n / a);
        else return complex<double>(0, 0);
    }

    /** banked received signal z_{ell,n} from the paper*/
    inline complex<double> z(int ell, int n) const {
        return zn(b * n + ell);
    }

    /** The sequence g_{ell,n} from the paper */
    inline complex<double> gfunc(int ell, int n) const {
        return conj(this->g.pulse(-(n - 1) * this->Delta - this->taumin + (ell * this->Ts) / a));
    }

    /** The sequence bk from the paper */
    virtual complex<double> v(int k) const {
        complex<double> sum(0, 0);
        for (unsigned int ell = 0; ell < b; ell++) sum += h(ell, k);
        //cout << b << ", " << sum << endl;
        return sum*this->Ts; //rescale discrete inner products by sample period.  This is not in the original paper, but it results in the same estimator
    }

    /** The h_{\ell,k} sequence from the paper, computed by convolution*/
    virtual complex<double> h(int ell, int k) const {
        double A = 1 - (this->g.tmax() + this->taumin) / this->Delta + ((double) ell) / b;
        double B = 1 - (this->g.tmin() + this->taumin) / this->Delta + ((double) ell) / b;
        int Bprime = (int) ceil((k - B + ell * n0) / a);
        int Aprime = (int) floor((k - A + ell * n0) / a);
        complex<double> sum(0, 0);
        int Bpp = a * Bprime - ell * n0;
        int App = ((int) a) * Aprime - ell * n0;
        for (int n = Bpp; n <= App; n += a) sum += z(ell, n) * gfunc(ell, k - n);
        return sum;
    }

    /** fill vector Zgrid with discretized Z function */
    virtual void fillZgrid() = 0;

    /** fill vector Ygridmag with discretized magnitude of the Y function */
    virtual void fillYgrid() = 0;

    //output for testing

    vector<double> getZgrid() {
        return Zgrid;
    }

    vector<double> getYgrid() {
        return Ygridmag;
    }

    vector<double> getSSgrid() {
        return SSgrid;
    }

protected:
    const unsigned int p;
    const unsigned int q;
    unsigned int a;
    unsigned int b;
    int n0;
    int m0;
    /** Size of the search grid */
    const unsigned int K;
    vector<double> Zgrid;
    vector<double> Ygridmag;
    vector<double> SSgrid;

};

/**
 * Fast polyphase computer for bk sequence.  Stores the bk sequence for fast access.
 */
template <class Pulse>
class FastAlgorithm : public FastAlgorithmNoStorev<Pulse> {
public:

    FastAlgorithm(const vector<int>& P,
            const vector<int>& D,
            const vector<complex<double> >& pilots,
            const Pulse& g,
            const double T,
            const double Ts,
            const double taumin,
            const double taumax,
            const unsigned int c,
            const unsigned int p,
            const unsigned int q,
            const double brenttol = DEFAULTBRENTTOL) :
    FastAlgorithmNoStorev<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q, brenttol) {
        int Dmin = *std::min_element(D.begin(), D.end());
        int Dmax = *std::max_element(D.begin(), D.end());
        int Pmin = *std::min_element(P.begin(), P.end());
        int Pmax = *std::max_element(P.begin(), P.end());
        vstore.resize(this->K + c*(max(Dmax, Pmax) - min(Dmin, Pmin)), 0.0);
        boffset = 1 + c*(min(Dmin, Pmin) + 1);
    };

    virtual double estimate(const vector<complex<double>>& r) {
        this->setupr(r);
        for (unsigned int k = 0; k < vstore.size(); k++) vstore[k] = FastAlgorithmNoStorev<Pulse>::v(k + boffset);
        double tautilde = this->coarseMaximiseSS();
        double tauhat = this->refineCoarseEstimate(tautilde);
        return tauhat;
    }

    //b takes from array
    virtual complex<double> v(int k) const {
        //if(k - boffset < 0 || k - boffset >= bstore.size()) std::cout << "Indexing out of bounds!" << std::endl;
        return vstore[k - boffset];
    }

protected:
    vector<complex<double> > vstore;
    int boffset;

};

template <class Pulse>
class Direct : public FastAlgorithm<Pulse> {
public:

    Direct(const vector<int>& P,
            const vector<int>& D,
            const vector<complex<double>>& pilots,
            const Pulse& g,
            const double T,
            const double Ts,
            const double taumin,
            const double taumax,
            const unsigned int c,
            const unsigned int p,
            const unsigned int q,
            const double brenttol = DEFAULTBRENTTOL) :
    FastAlgorithm<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q, brenttol) {
    }

    /** Return the values of Z computed on the grid taumin to taumax by Ts/c.  Direct computation. */
    virtual void fillZgrid() {
        for (unsigned int k = 1; k <= this->K; k++) {
            double sum = 0.0;
            for (int i : this->D)
                sum += std::abs(this->v(k + this->c * (i + 1)));
            this->Zgrid[k - 1] = sum;
        }
    }

    /* Return the values of Y computed on the grid taumin to taumax by T/c. Direct convolution */
    virtual void fillYgrid() {
        for (unsigned int k = 1; k <= this->K; k++) {
            complex<double> sum = 0.0;
            for (unsigned int i = 0; i < this->P.size(); i++)
                sum += this->v(k + this->c * (this->P[i] + 1)) * conj(this->pilots[i]);
            this->Ygridmag[k - 1] = std::abs(sum);
        }
    }

    virtual string name() const {
        return "Direct";
    }
};

/**
 * Recursive estimator of time offset.  The data symbol indices must now be 
 * contiguous.  Exception is thrown if they are not.
 */
template <class Pulse>
class Recursive : public Direct<Pulse> {
public:
    Recursive(const vector<int>& P,
            const vector<int>& D,
            const vector<complex<double>>& pilots,
            const Pulse& g,
            const double T,
            const double Ts,
            const double taumin,
            const double taumax,
            const unsigned int c,
            const unsigned int p,
            const unsigned int q,
            const double brenttol = DEFAULTBRENTTOL) :
    Direct<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q, brenttol),
    Dmin(*std::min_element(D.begin(), D.end())),
    Dmax(*std::max_element(D.begin(), D.end())) {
        //check D is contiguous and sorted
        for (unsigned int i = 0; i < D.size() - 1; i++) {
            if (D[i] != (D[i + 1] - 1)) {
                std::ostringstream str;
                str << "Exception constructing Recursive. The data symbols must be contiguous." << endl;
                throw str.str();
            }
        }
    }

    /* 
     * Return the values of Z computed on the grid taumin to taumax by T/c. 
     * Uses recursive algorithm that only applies when D is contiguous.
     */
    virtual void fillZgrid() {
        for (unsigned int k = 1; k <= this->c; k++) {
            double sum = 0.0;
            for (int i : this->D) sum += abs(this->v(k + this->c * (i + 1)));
            this->Zgrid[k - 1] = sum;
            for (int m = 0; k + (m + 1) * this->c <= this->K; m++)
                this->Zgrid[k - 1 + (m + 1) * this->c] = this->Zgrid[k - 1 + m * this->c] - abs(this->v(k + m * this->c + (Dmin + 1) * this->c)) + abs(this->v(k + m * this->c + (Dmax + 1) * this->c + this->c));
        }
    }

    virtual string name() const {
        return "Recursive";
    }
    
    protected:
        /** Minimum data symbol index */
        const int Dmin;
        /** Maximum data symbol index */
        const int Dmax;
};

/**
 * Estimator of time offset that only uses pilot symbols
 */
template <class Pulse>
class PilotsOnly : public Direct<Pulse> {
public:
    PilotsOnly(const vector<int>& P,
            const vector<complex<double>>& pilots,
            const Pulse& g,
            const double T,
            const double Ts,
            const double taumin,
            const double taumax,
            const unsigned int c,
            const unsigned int p,
            const unsigned int q,
            const double brenttol = DEFAULTBRENTTOL) :
    Direct<Pulse>(P, P, pilots, g, T, Ts, taumin, taumax, c, p, q, brenttol) {
    }

    /** Does nothing! There are only pilots */
    inline virtual void fillZgrid() { }
    
    /** Return zero! There are only pilots */
    virtual double Z(double tau) const {
        return 0;
    }
    
    /** 
     * This replaces the pilots, copy new ones.  This is really just a hack to make simulations
     * easier!  Ideally this would not exist and pilots would be const
     */
    virtual void setPilots(const vector<complex<double>>& pilots){
        if(this->pilots.size() != pilots.size()) throw "You shouldn't be changing the number of pilots";
        for(int i = 0; i < pilots.size(); i++) this->pilots[i] = pilots[i];
    }
    
    virtual string name() const {
        return "PilotsOnly";
    }
   
};

#endif	/* FASTALGORITHM_H */

