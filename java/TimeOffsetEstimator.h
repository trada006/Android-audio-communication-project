/* 
 * File:   TimeOffsetEstimator.h
 * Author: Robby McKilliam
 *
 * Created on 23 January 2013, 10:58 AM
 */

#ifndef TIMEOFFSETESTIMATOR_H
#define	TIMEOFFSETESTIMATOR_H

#include "FinitePulse.h"
#include <complex>
#include <vector>
#include <algorithm>

#define DEFAULTBRENTTOL 1e-7

using namespace std;

/** Interface for time offset estimators  */
class TimeOffsetEstimator {
public:
    virtual ~TimeOffsetEstimator() {
    };
    /** Run the estimator. Return time offset estimate */
    virtual double estimate(const vector<complex<double>>& r) = 0;
    /** Return a string with the name of this estimator */
    virtual string name() const = 0;
 };

/** Naive implementation of the time offset estimator that computes the objective function directly */
template <class Pulse> //P should be a FinitePulse
class Naive : public TimeOffsetEstimator {
public:

    Naive(const vector<int>& P_,
            const vector<int>& D_,
            const vector<complex<double>>& pilots_,
            const Pulse& g_,
            const double T_,
            const double Ts_,
            const double taumin_,
            const double taumax_,
            const unsigned int c_,
            const double brenttol_ = DEFAULTBRENTTOL) :
    L(D_.size() + P_.size()),
    Delta(T_ / c_),
    P(P_),
    D(D_),
    pilots(pilots_),
    g(g_),
    T(T_),
    Ts(Ts_),
    taumin(taumin_),
    taumax(taumax_),
    c(c_),
    brenttol(brenttol_*T){
        if(taumax <= taumin) throw "Maximum time offset taumax must be larger than minimum time offset taumin";
    };

    virtual double estimate(const vector<complex<double>>& r) {
        setupr(r);
        double tautilde = coarseMaximiseSS();
        double tauhat = refineCoarseEstimate(tautilde);
        return tauhat;
    }

    /** Obtain a coarse estimate of the time offset */
    virtual double coarseMaximiseSS() {
        double maxSS = -1.0;
        double taubest = taumin;
        for (double tau = taumin; tau <= taumax; tau += Delta) {
            double thisSS = SS(tau);
            if (maxSS < thisSS) {
                maxSS = thisSS;
                taubest = tau;
            }
        }
        return taubest;
    }

    /** Refines the coarse estimate. */
    virtual double refineCoarseEstimate(double tautilde) const {
        double a = tautilde - Delta;
        double c = tautilde + Delta;
        auto f = [this] (double tau){return -this->SS(tau);}; //function to minimise
        Brent opt(f, a, tautilde, c, brenttol);
        return opt.xmin();
    }

    /** Inner product between the received signal and pulse g */
    virtual complex<double> m(double tau) const {
        int A = (int) ceil((g.tmin() + tau) / Ts);
        int B = (int) floor((g.tmax() + tau) / Ts);
        complex<double> sum(0, 0);
        for (int n = A; n <= B; n++) sum += r(n) * conj(g.pulse(n * Ts - tau));
        return sum*Ts; //rescale discrete inner products by sample period.  This is not in the original paper, but it results in the same estimator
    }

    /** The Y function computing a correlation with the data and pilots */
    virtual complex<double> Y(double tau) const {
        complex<double> sum(0, 0);
        for (unsigned int i = 0; i < P.size(); i++) sum += m((P[i] + 1) * T + tau) * conj(pilots[i]);
        return sum;
    }

    /** The amplitude accumulating Z function */
    virtual double Z(double tau) const {
        double sum = 0.0;
        for (int i : D) {
            sum += abs(m((i + 1) * T + tau));
        }
        return sum;
    }

    /** The objective function */
    virtual inline double SS(double tau) const {
        return Z(tau) + abs(Y(tau));
    }

    virtual string name() const {
        return "Naive";
    }

    const Pulse g;
    const double T;
    const double Ts;
    const double taumin;
    const double taumax;
    const unsigned int c;
    const double brenttol;
    /** Total number of transmitted symbols */
    const unsigned int L;
    /** Grid search width */
    const double Delta;
    
protected:

    /** Pointer to current received signal */
    const complex<double>* rmem;
    /** size of current received signal */
    unsigned int rmemsize;

    /** The received sequence, zeros returned for unknown values */
    inline complex<double> r(int n) const {
        if (n > 0 && n <= rmemsize) return rmem[n - 1];
        else return complex<double>(0, 0);
    }
    
    virtual void setupr(const vector<complex<double>>& r) {
        rmem = r.data();
        rmemsize = r.size();
    }

    const vector<int> P;
    const vector<int> D;
    vector<complex<double>> pilots;
    
};

#endif	/* TIMEOFFSETESTIMATOR_H */

