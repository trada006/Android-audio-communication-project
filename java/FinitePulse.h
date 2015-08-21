/* 
 * File:   FinitePulse.h
 * Author: Robby McKilliam
 *
 * Created on 23 January 2013, 1:10 PM
 */

#ifndef FINITEPULSE_H
#define	FINITEPULSE_H

#include <complex>
#include <vector>
#include <iostream>
#include "Util.h"

using namespace std;

class FinitePulse {
public:

    virtual ~FinitePulse() {
    };
    /** pulse is zero for t less than tmin */
    virtual double tmin() const = 0;
    /** pulse is zero for t more than tmax */
    virtual double tmax() const = 0;
    /** pulse period */
    virtual double T() const = 0;
    /** The transmit pulse */
    virtual complex<double> pulse(double t) const = 0;
};

/** Truncated sinc pulse goes out to the numzeros zero in both positive and negative directions */
class TruncatedSincPulse : public FinitePulse {
public:

    TruncatedSincPulse(double T, unsigned int numzeros) : tmin_(-T*numzeros), tmax_(T*numzeros), T_(T), sqrtT(sqrt(T)) {
    };
    TruncatedSincPulse() = delete; //no default constructor

    /** The truncated sinc pulse */
    virtual complex<double> pulse(double t) const {
        if (t > tmin_ && t < tmax_) return sinc(t / T_) / sqrtT;
        else return 0.0;
    }

    virtual double tmin() const {
        return tmin_;
    }

    virtual double tmax() const {
        return tmax_;
    }

    virtual double T() const {
        return T_;
    }

protected:
    const double tmin_;
    const double tmax_;
    const double T_;
    const double sqrtT;
};

/** 
 * Takes a pulse of finite duration and normalises it to have energy 1.
 * Uses numerical integration to compute the normalising constant
 */
template <class P> //P should be a FinitePulse
class NormalisedFinitePulse : public FinitePulse {
public:

    NormalisedFinitePulse(const P& p_) : p(p_), tmin_(p_.tmin()), tmax_(p_.tmax()), T_(p_.T()) {
        auto f = [this](double t) -> double{return std::norm(this->p.pulse(t));};
        double energy = trapezoidal(f, tmin_, tmax_, 100000); //last number is integration steps
        normalisingconstant = sqrt(energy);
    };
    //~NormalisedFinitePulse() { delete p; }
    //NormalisedFinitePulse( const NormalisedFinitePulse& other ) = delete; //no copy
    //NormalisedFinitePulse& operator =(const NormalisedFinitePulse&) = delete; //no assignment
    //NormalisedFinitePulse() = delete; //no default constructor

    virtual double tmin() const {
        return tmin_;
    }

    virtual double tmax() const {
        return tmax_;
    }

    virtual double T() const {
        return T_;
    }

    /** The normalised pulse */
    virtual complex<double> pulse(double t) const {
        return p.pulse(t) / normalisingconstant;
    }

protected:
    const P p;
    const double tmin_;
    const double tmax_;
    const double T_;
    double normalisingconstant;

};

class TruncatedRootRaisedCosine : public FinitePulse {
public:
    TruncatedRootRaisedCosine(double T, double beta, double duration, bool dodgyflag) : T_(T), beta(beta)  {
        const double stepsize = 0.01; //step taken whilst looking for zeros (this will work only if zeros are atleast stepsize apart)
        double c = 0.0;
        int csign = 1;
        while(c < duration/T/2) {
            while (signum(rootraisedcosine(c)) == csign) c += stepsize;
            csign = -csign;
        }
        double nt = Bisection([this] (double x) { return rootraisedcosine(x); }, c - stepsize, c, 1e-7, 100).zero();
        tmin_ = -nt*T;
        tmax_ = nt*T;
        double dur = tmax_ - tmin_;
        //std::cout << dur << " but request duration is " << duration << std::endl;
        if( dur > 2*duration || dur < duration/2 ) { //check that the duration obtained is reasonable
            std::ostringstream str;
            str << "Something when wrong with the bisection method, duration is ";
            str << dur << " but request duration is " << duration;
            throw str.str();
        }
    }
    TruncatedRootRaisedCosine(double T, double beta, int numzeros) :
    T_(T), beta(beta) {
        const double stepsize = 0.01; //step taken whilst looking for zeros (this will work only if zeros are atleast stepsize apart)
        double c = 0.0;
        int csign = 1;
        for (int i = 1; i <= numzeros; i++) {
            while (signum(rootraisedcosine(c)) == csign) c += stepsize;
            csign = -csign;
        }

        double nt = Bisection([this] (double x) { return rootraisedcosine(x);}, c - stepsize, c, 1e-7, 100).zero();
        tmin_ = -nt*T;
        tmax_ = nt*T;
    }
    TruncatedRootRaisedCosine() = delete; //no default constructor

    virtual complex<double> pulse(double t) const {
        if (t > tmin_ && t < tmax_) return rootraisedcosine(t / T_);
        else return 0.0;
    }

    virtual double tmin() const {
        return tmin_;
    }

    virtual double tmax() const {
        return tmax_;
    }

    virtual double T() const {
        return T_;
    }

    /** A root raised cosine with period 1 and rolloff beta */
    double rootraisedcosine(double t) const {
        double abst = fabs(t); //pulse is symmetric about zero, so just use magnitude of t
        if (abst < 5e-3) { //second order expansion if t is near zero
            double term0 = 1 + beta * (4 / pi - 1);
            double term2 = (cub((beta - 1) * pi) + 96 * beta * beta * (4 * beta + pi - beta * pi) - 12 * beta * sqr(pi + beta * pi)) / 6 / pi;
            return term0 + term2 * abst*abst;
        }
        if (fabs(abst - 1.0 / 4 / beta) < 5e-4) { //first order expansion if t is near 1/(4beta)
            double a = (1 + beta) * pi / 4 / beta;
            double term0 = beta * (sin(a) - 2 * cos(a) / pi);
            double term1 = beta * ((12 * beta + pi * pi) * cos(a) + 2 * (1 - 2 * beta) * pi * sin(a)) / pi;
            return term0 + term1 * (abst - 1.0 / 4 / beta);
        } else { //otherwise use direct formula
            double a = pi*t;
            double b = beta*t;
            return (sin(a * (1 - beta)) + 4 * b * cos(a * (1 + beta))) / (a * (1 - 16 * b * b));
        }
    }

protected:
    double tmin_;
    double tmax_;
    const double T_;
    const double beta;
};

#endif	/* FINITEPULSE_H */

