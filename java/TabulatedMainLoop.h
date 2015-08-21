/* 
 * File:   TabulatedMainLoop.h
 * Author: Robby McKilliam
 *
 * Created on 4 March 2013, 11:38 PM
 */

#ifndef TABULATEDMAINLOOP_H
#define	TABULATEDMAINLOOP_H

#include <complex>

using namespace std;

/** Runs the main filter loop, this is where most time is spent */
static inline complex<double> mainFilterLoop(int mfrom, int mto, int mstep, int ifrom, int istep, const complex<double>* rmem, const complex<double>* pulsetable) {
    complex<double> sum(0, 0);
    //mainloop, this will greatly benefit from a mac instruction
    for (int m = mfrom, i = ifrom; m <= mto; m += mstep, i += istep)
        sum += rmem[m] * pulsetable[i];
    return sum;
}

/** 
 * Runs the main filter loop. This assumes you have arranged memory so that the pulsetable 
 * can be stepped through with unit increments.
 */
static inline complex<double> mainBankedFilterLoop(int mfrom, int mto, int mstep, const complex<double>* rmem, const complex<double>* pulsetable) {
    complex<double> sum(0, 0);
    //mainloop, this will greatly benefit from a mac instruction
    for (int m = mfrom, i = 0; m <= mto; m += mstep, i++)
        sum += rmem[m] * pulsetable[i];
    return sum;
}

#endif	/* TABULATEDMAINLOOP_H */

