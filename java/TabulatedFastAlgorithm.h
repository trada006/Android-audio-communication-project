/* 
 * File:   TabulatedFastAlgorithm.h
 * Author: Robby McKilliam
 *
 * Created on 20 February 2013, 6:36 PM
 */

#ifndef TABULATEDFASTALGORITHM_H
#define	TABULATEDFASTALGORITHM_H

#include <complex>
#include <vector>
#include <algorithm>
#include "TimeOffsetEstimator.h"
#include "FinitePulse.h"
#include "TabulatedMainLoop.h"
#include "FastAlgorithm.h"

#define DEFAULTTABLESIZE 100000

/**
 * Version of the polyphase TimeOffset estimator that tabulates the transmit pulse in a way to
 * allow fast access, but also fast, well pipeliped loops for filters.
 */
template <class Pulse>
class TabulatedDirect : public FastAlgorithm<Pulse> {
public:

    TabulatedDirect(const vector<int>& P,
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
            const unsigned int mintabs = DEFAULTTABLESIZE,
            const double brenttol = DEFAULTBRENTTOL) :
    FastAlgorithm<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q, brenttol) {
        //int d = lcm(q,c);
        int d = c * this->b; //we require the number of table elements per symbol to be a multiple of d
        tabls = (mintabs + 1) / d;
        int tablesizepersymbol = tabls*d;
        stepwidth = T / tablesizepersymbol;
        for (double t = g.tmin(); t <= g.tmax(); t += stepwidth) pulsetable.push_back(conj(g.pulse(t)));
    };

    /** Inner product between the received signal and pulse g */
    inline virtual complex<double> m(double tau) const {

        //step loop counters.
        int A = (int) ceil((this->g.tmin() + tau) / this->Ts);
        int B = (int) floor((this->g.tmax() + tau) / this->Ts);
        int nfrom = max(1, A);
        int nto = min((int) this->rmemsize, B);
        double startt = nfrom * this->Ts - tau;
        int ifrom = (int) round((startt - this->g.tmin()) / this->stepwidth);
        int istep = this->a*tabls;
        int mfrom = nfrom - 1;
        int mto = nto - 1;
        int mstep = 1;

        return mainFilterLoop(mfrom, mto, mstep, ifrom, istep, this->rmem, &pulsetable[0]) * this->Ts; //rescale discrete inner products by sample period.  This is not in the original paper, but it results in the same estimator
    }

    /** The h_{\ell,k} sequence from the paper, computed by convolution*/
    inline virtual complex<double> h(int ell, int k) const {

        //step loop counters.
        double A = 1 - (this->g.tmax() + this->taumin) / this->Delta + ((double) ell) / this->b;
        double B = 1 - (this->g.tmin() + this->taumin) / this->Delta + ((double) ell) / this->b;
        int Bprime = (int) ceil((k - B + ell * this->n0) / this->a);
        int Aprime = (int) floor((k - A + ell * this->n0) / this->a);
        int Bpp = this->a * Bprime - ell * this->n0;
        int App = this->a * Aprime - ell * this->n0;
        int nfrom = max((int) (ell / this->b + 1), Bpp);
        //int nto = min((int) (a * (rmemsize + 1) / b), App);
        int nto = App;
        double startt = -(k - nfrom - 1) * this->Delta - this->taumin + (ell * this->Ts) / this->a;
        int ifrom = (int) round((startt - this->g.tmin()) / stepwidth);
        int istep = this->a*this->b*tabls;
        int mfrom = (this->b * nfrom + ell) / this->a - 1;
        //int mto = (b * nto + ell) / a - 1;
        int mto = min((int)this->rmemsize,(int)((this->b*nto+ell)/this->a))-1;
        int mstep = this->b;

        return mainFilterLoop(mfrom, mto, mstep, ifrom, istep, this->rmem, &pulsetable[0]);
    }

    /** Return the values of Z computed on the grid taumin to taumax by Ts/c.  Direct computation. */
    virtual void fillZgrid() {
        for (unsigned int k = 1; k <= this->K; k++) {
            double sum = 0.0;
            for (int i = 0; i < this->D.size(); i++) sum += std::abs(this->v(k + this->c * (this->D[i] + 1)));
            //for (int i : this->D) sum += std::abs(bfunc(k + this->c * (i + 1)));
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
        return "TabulatedDirect";
    }

protected:
    /** Table stores values of the transmit pulse g.pulse */
    std::vector<complex<double>> pulsetable;
    /**  Width in t between elements in the table. */
    double stepwidth;
    /** Multiplier for indexing the table */
    int tabls;

};

/**
 * Recursive estimator of time offset.  The data symbol indices must now be 
 * contiguous.  Exception is thrown if they are not.
 */
template <class Pulse>
class TabulatedRecursive : public TabulatedDirect<Pulse> {
public:

    TabulatedRecursive(const vector<int>& P,
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
            const unsigned int mintabs = DEFAULTTABLESIZE,
            const double brenttol = DEFAULTBRENTTOL) :
    TabulatedDirect<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q, mintabs, brenttol),
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
            for (int i = 0; i < this->D.size(); i++) sum += abs(this->v(k + this->c * (this->D[i] + 1)));
            this->Zgrid[k - 1] = sum;
            for (int m = 0; k + (m + 1) * this->c <= this->K; m++)
                this->Zgrid[k - 1 + (m + 1) * this->c] = this->Zgrid[k - 1 + m * this->c] - abs(this->v(k + m * this->c + (Dmin + 1) * this->c)) + abs(this->v(k + m * this->c + (Dmax + 1) * this->c + this->c));
        }
    }

    virtual string name() const {
        return "TabulatedRecursive";
    }

protected:
    /** Minimum data symbol index */
    const int Dmin;
    /** Maximum data symbol index */
    const int Dmax;
};

/**
 * Recursive estimator of time offset.  Does no refinement, i.e. no Brent's method.  Just here
 * for benchmarking.
 */
template <class Pulse>
class TabulatedRecursiveNoRefine : public TabulatedRecursive<Pulse> {
public:

    TabulatedRecursiveNoRefine(const vector<int>& P,
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
            const unsigned int mintabs = DEFAULTTABLESIZE) :
    TabulatedRecursive<Pulse>(P, D, pilots, g, T, Ts, taumin, taumax, c, p, q, mintabs) {
    }

    /** No refinement. */
    virtual double refineCoarseEstimate(double tautilde) const {
        return tautilde;
    }
    
    virtual string name() const {
        return "TabulatedRecursiveNoRefine";
    }

};

/**
 * Recursive estimator of time offset.  The data symbol indices must now be 
 * contiguous.  Exception is thrown if they are not.
 */
template <class Pulse>
class TabulatedPilotsOnly : public TabulatedDirect<Pulse> {
public:
    TabulatedPilotsOnly(const vector<int>& P,
            const vector<complex<double> >& pilots,
            const Pulse& g,
            const double T,
            const double Ts,
            const double taumin,
            const double taumax,
            const unsigned int c,
            const unsigned int p,
            const unsigned int q,
            const unsigned int mintabs = DEFAULTTABLESIZE,
            const double brenttol = DEFAULTBRENTTOL) :
    TabulatedDirect<Pulse>(P, P, pilots, g, T, Ts, taumin, taumax, c, p, q, mintabs, brenttol) {
    }

    /** Does nothing! There are only pilots */
    inline virtual void fillZgrid() { }
    
    /** Return zero! There are only pilots */
    inline virtual double Z(double tau) const {
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
        return "TabulatedPilotsOnly";
    }
   
};

#endif	/* TABULATEDFASTALGORITHM_H */

