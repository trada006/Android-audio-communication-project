/* 
 * File:   QAMTimeOffsetEstimator.h
 * Author: Robby McKilliam
 *
 * Created on 27 January 2015, 9:58 AM
 */

#ifndef QAMTIMEOFFSETESTIMATOR_H
#define	QAMTIMEOFFSETESTIMATOR_H

#include "Util.h"
#include "TabulatedFastAlgorithm.h"

/**
 * Modified timeoffset estimator for QAM symbols (i.e. varying amplitude. The objective function
 * looks a bit like the subspace detectors of Scharf, except with the nonlinear time offset 
 * component. One advantage of this objective function is that it remains differentiable
 * everywhere. This might make a theoretical analysis more tractable.
 */
template <class Pulse>
class QAMDirect : public TabulatedDirect<Pulse> {
public:

    QAMDirect(const vector<int>& P,
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
    Ep(std::real(inner_product(pilots,pilots))) {}
    
    virtual string name() const {
        return "QAMDirect";
    }
    
    /** Energy of the pilot symbols */
    const double Ep;

    /** The Y function computing a correlation with the data and pilots */
    virtual complex<double> Y(double tau) const {
        complex<double> sum(0, 0);
        for (unsigned int i = 0; i < this->P.size(); i++) sum += this->m((this->P[i] + 1) * this->T + tau) * conj(this->pilots[i]);
        return sum;
    }

    /** The amplitude accumulating Z function */
    virtual double Z(double tau) const {
        double sum = 0.0;
        for (int i : this->D) {
            sum += norm(this->m((i + 1) * this->T + tau));
        }
        return sum;
    }
    
    /** Return the values of Z computed on the grid taumin to taumax by Ts/c.  Direct computation. */
    virtual void fillZgrid() {
        for (unsigned int k = 1; k <= this->K; k++) {
            double sum = 0.0;
            for (int i = 0; i < this->D.size(); i++) sum += std::norm(this->v(k + this->c * (this->D[i] + 1)));
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
            this->Ygridmag[k - 1] = std::norm(sum)/Ep;
        }
    }

    /** The objective function */
    virtual inline double SS(double tau) const {
        return this->Z(tau) + norm(this->Y(tau))/Ep;
    }
    
};

#endif	/* QAMTIMEOFFSETESTIMATOR_H */

