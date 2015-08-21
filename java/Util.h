/* 
 * File:   Util.h
 * Author: Robby McKilliam
 *
 * Created on 23 January 2013, 2:55 PM
 */

#ifndef UTIL_H
#define	UTIL_H

#include <complex>
#include <functional>
#include <vector>
#include <algorithm> 

static constexpr double pi = 3.141592653589793238463;

/** The sign of x, zero if x is zero. */
inline static double signum(double x) {
    return (x > 0) - (x < 0);
}

/** Square of x */
inline static double sqr(double x) {
    return x*x;
}

/** Cube of x */
inline static double cub(double x) {
    return x * x*x;
}

/** sinc function */
static double sinc(double t) {
    if (fabs(t) < 5e-3) return 1.0 - t * t * (1.0 / 6 - 1.0 / 120 * t * t);
    else return sin(pi * t) / (pi * t);
}

/** Trapezoidal integration of a function f mapping double to double */
static double trapezoidal(std::function<double(double)> f, double a, double b, unsigned int N) {
    double del = (b - a) / N;
    double inner = 0.0;
    for (unsigned int n = 1; n <= N - 1; n++) inner += 2 * f(a + n * del);
    return del / 2 * (inner + f(a) + f(b));
}

/** Greatest common divisor between two integers */
static int gcd(int a, int b) {
    if (b == 0) return abs(a);
    else return gcd(abs(b), abs(a) % abs(b));
}

template <class T>
static std::vector<T> concatenate(const std::vector<T>& A, const std::vector<T>& B) {
    std::vector<T> C;
    C.reserve(A.size() + B.size());
    for (T a : A) C.push_back(a);
    for (T b : B) C.push_back(b);
    return C;
}

/** 
 * The Extended Euclidean algorithm applied to two integers. Returns integers x, y and gcd 
 * such that ax + by = gcd.
 */
static void extended_gcd(int a, int b, int& gcd, int& x, int& y) {
    x = 0, y = 1;
    int u = 1, v = 0, m, n, q, r;
    gcd = b;
    while (a != 0) {
        q = gcd / a;
        r = gcd % a;
        m = x - u*q;
        n = y - v*q;
        gcd = a;
        a = r;
        x = u;
        y = v;
        u = m;
        v = n;
    }
}

/** The Hermitian inner product between two complex vectors. */
static std::complex<double> inner_product(const std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& y){
    int N = x.size();
    if(N != y.size()) throw "Vectors must be the same size for inner product";
    std::complex<double> sum(0);
    for(int i = 0; i < N; i++){
        sum += x[i]*std::conj(y[i]);
    }
    return sum;
}

/**
 *Performs a 1-dimensional minimization.
 *It implements Brent's method which combines a golden-section
 *search and parabolic interpolation. 
 *
 *@param f           The function to minimise 
 *@param  ax         Left endpoint of initial interval
 *@param  cx         Right endpoint of initial interval
 *@param  bx         bx must satisfy ax < bx < cx and f(ax) > f(bx) < f(cx)
 *@param  tol         Desired length of the interval in which the minimum will be determined to lie (default 100)
 *@params ITRMAX  maximum number of iterations before terminating (default 1e-8)
 */
class Brent {
public:

    Brent(std::function<double(double) > f, double ax, double bx, double cx, double tol = 1e-8, unsigned int ITRMAX = 100) {
        double e = 0.0;
        double d = 0.0;
        double a = (ax < cx ? ax : cx);
        double b = (ax > cx ? ax : cx);
        x = bx;
        double w = bx;
        double v = bx;
        double fw = f(x);
        double fv = fw;
        fx = fw;
        for (unsigned int iter = 0; iter < ITRMAX; iter++) {
            const double xm = 0.5 * (a + b);
            const double tol1 = tol * fabs(x) + ZEPS;
            const double tol2 = 2.0 * tol1;
            if (fabs(x - xm) <= tol2 - 0.5 * (b - a)) {
                return; //(fx, x)
            }
            if (fabs(e) > tol1) {
                const double r = (x - w)*(fx - fv);
                double q = (x - v)*(fx - fv);
                double p = (x - v) * q - (x - v) * r;
                q = 2.0 * (q - r);
                if (q > 0.0) p = -p;
                q = fabs(q);
                const double etemp = e;
                e = d;
                if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) { //golden step
                    e = (x >= xm) ? a - x : b - x;
                    d = C*e;
                } else { //parabolic step
                    d = p / q;
                    const double u = x + d;
                    if (u - a < tol2 || b - u < tol2) d = sign(tol1, xm - x);
                }
            } else {
                e = (x >= xm) ? a - x : b - x;
                d = C*e;
            }
            const double u = (fabs(d) >= tol1 ? x + d : x + sign(tol1, d));
            const double fu = f(u);
            if (fu <= fx) {
                if (u >= x) a = x;
                else b = x;
                v = w;
                w = x;
                x = u;
                fv = fw;
                fw = fx;
                fx = fu;
            } else {
                if (u < x) a = u;
                else b = u;
                if (fu <= fw || w == x) {
                    v = w;
                    w = u;
                    fv = fw;
                    fw = fu;
                } else if (fu <= fv || v == x || v == w) {
                    v = u;
                    fv = fu;
                }
            }
        }
        std::cout << "Warning: Brent's method reached the maximum number " << ITRMAX << " iterations" << std::endl;
        //return (fx, x)
    }

    /** The value of the function at the minimum */
    double fmin() const {
        return fx;
    }

    /** The minimiser of the function */
    double xmin() const {
        return x;
    }

    inline double sign(double a, double b) const {
        return fabs(a) * signum(b);
    }

private:
    static constexpr double ZEPS = 1e-10; //close to machine precision
    static constexpr double C = (3.0 - sqrt(5.0)) / 2.0;
    double fx;
    double x;

};

/** 
 *Search for a zero of the function f(x) in the interval [a, b].
 *Finds a solution such that |x0 - x| < tol where f(x0) = 0.
 *Uses the bisection method, only guaranteed to converge if there is a unique zero 
 *between a and b and f is continuous and sign(f(a)) = - sign(f(b))
 *
 *@param f           The function to zero 
 *@param  a         Left endpoint of initial interval
 *@param  c         Right endpoint of initial interval
 *@param  b         bx must satisfy a < b < c and f(a) > f(b) < f(c)
 *@param  tol         Desired length of the interval in which the minimum will be determined to lie (default 1e-6)
 *@params ITRMAX  maximum number of iterations before terminating (default 100)
 */
class Bisection {
public:

    Bisection(std::function<double(double) > f, double ax, double bx, double tol = 1e-8, unsigned int ITRMAX = 100) {
        double a = ax;
        double b = bx;
        for (unsigned int i = 1; i <= ITRMAX; i++) {
            double c = (a + b) / 2;
            double fc = f(c);
            if (fc == 0 || fabs(a - b) / 2 < tol) {
                xzero = c;
                return;
            }
            if (signum(fc) == signum(f(a))) a = c;
            else b = c;
        }
        std::cout << "Warning: Bisection failed. Reached the maximum number " << ITRMAX << " iterations" << std::endl;
    }

    /** The value of x where f(x) = 0 */
    double zero() const {
        return xzero;
    }

protected:
    double xzero;

};

/** Generates an m-sequence of length 2^N-1. */
static std::vector<int> msequence(int N){
    if(N < 3 || N > 14) throw "Only m-sequence with N > 2 and N < 10 and available.";
    
    //initial taps
    std::vector<int> taps;
    if(N == 3) {
        int array[3] = {0,1,1};
        for(int i = 0; i < N; i++) taps.push_back(array[i]);
    }
    else if(N == 4){
        int array[4] = {0,0,1,1};
        for(int i = 0; i < N; i++) taps.push_back(array[i]);
    }
    else if(N == 5){
        int array[5] = {0,0,1,0,1};
        for(int i = 0; i < N; i++) taps.push_back(array[i]);
    }
    else if(N == 6){
        int array[6] = {0,0,0,0,1,1};
        for(int i = 0; i < N; i++) taps.push_back(array[i]);
    }
    else if(N == 7){
        int array[7] = {0,0,0,1,0,0,1};
        for(int i = 0; i < N; i++) taps.push_back(array[i]);
    }
    else if(N == 8){
        int array[8] = {0,0,0,1,1,1,0,1};
        for(int i = 0; i < N; i++) taps.push_back(array[i]);
    }
    else if(N == 9){
        int array[9] = {0,0,0,0,1,0,0,0,1};
        for(int i = 0; i < N; i++) taps.push_back(array[i]);
    }
    else if(N == 10){
      int array[10] = {0,0,0,0,0,0,1,0,0,1};
        for(int i = 0; i < N; i++) taps.push_back(array[i]);
    }
    else if(N == 11){
      int array[11] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1};
        for(int i = 0; i < N; i++) taps.push_back(array[i]);
    }
    else if(N == 12){
      int array[12] = {0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1};
        for(int i = 0; i < N; i++) taps.push_back(array[i]);
    }
    else if(N == 13){
      int array[13] = {0,0,0,0,0,0,0,0,1,1,0,1,1};
        for(int i = 0; i < N; i++) taps.push_back(array[i]);
    }
    else if(N == 14){
      int array[14] = {0,0,0,1,0,0,0,1,0,0,0,0,1,1};
        for(int i = 0; i < N; i++) taps.push_back(array[i]);
    }
    
    int M = (1<<N) - 1; //M = 2^N-1
    std::vector<int> m(N,1); //vector of 1's of length N
    std::vector<int> regout(M); //output m sequence

    for(int ind = 0; ind < M; ind++){
        int buf = 0;
        for(int i = 0; i < N; i++) buf += taps[i]*m[i];
        buf = buf % 2; //xor
        std::rotate(m.begin(), m.end()-1, m.end()); //shift
        m[0] = buf;
        regout[ind] = m[N-1];
    }
    
    return regout;
 
}

/**Generates a pilot sequence of length N^2 - 1 from an m-sequence by replacing 0's with 1's.*/
static std::vector<std::complex<double>> pilotmsequence(int N){
    std::vector<int> ms = msequence(N);
    std::vector<std::complex<double>> pilots(ms.size());
    for(int i = 0; i < ms.size(); i++){
        if(ms[i]==0) pilots[i] = -1;
        else pilots[i] = 1;
    }
    return pilots;
}

#endif	/* UTIL_H */

