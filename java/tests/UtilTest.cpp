/* 
 * File:   UtilTest.cpp
 * Author: Robby McKilliam
 *
 * Created on 24/01/2013, 11:21:20 AM
 */

#include <stdlib.h>
#include <iostream>
#include "Util.h"

double sqrint(double a, double b) { return b*b*b/3 - a*a*a/3; }
void trapezoidalSquareTest() {
    double tol = 1e-6;
    std::cout << "Trapezoidal Square Test ... ";
    unsigned int N = 5000;
    bool pass = abs(trapezoidal(sqr, -3.0, 4.0, N) - sqrint(-3.0,4.0)) < tol;   
    if(pass) std::cout << "PASS" << std::endl;
    else {
        std::cout << "FAIL" << std::endl;
        std::cout << "expected " << sqrint(-3.0,4.0) << ", but was " << trapezoidal(sqr, -3.0, 4.0, N) << std::endl;
    }
}

void gcdTest() {
    std::cout << "gcd ... ";
    bool pass = gcd(10,2)==2;
    pass &= gcd(3458,4864)==38;
    if(pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void extended_gcd_Test() {
    std::cout << "extended gcd ... ";
    int gcd, x, y;
    bool pass = true;
    extended_gcd(10,2,gcd,x,y);
    pass &= (gcd==2)&&((10*x + 2*y) == gcd);
    extended_gcd(3458,4864,gcd,x,y);
    pass &= (gcd==38)&&((3458*x + 4864*y) == gcd);
    if(pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void fzeroLinearTest() {
    std::cout << "Bisection linear ... ";
    double tol = 1e-7;
    auto f = [] (double x) { return x; };
    double x = Bisection(f, -11.0,8.0,tol).zero();
    bool pass = fabs(0.0 - f(x)) <  tol;
    pass &= fabs(0.0 - x) <  tol;
    if(pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
  }

void fzeroCubicTest() {
    std::cout << "Bisection cubic ... ";
    double tol = 1e-7;
    auto f = [] (double x) { return x*x*(x-1); };
    double x = Bisection(f, 0.5,1.7,tol).zero();
    bool pass = true;
    pass &= fabs(1.0 - x) <  tol;
    if(pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
  }

void msequenceTest() {
    std::cout << "m-sequence ... ";
    bool pass = true;
    std::vector<int> s = msequence(3);
    int exp3[7] = {1,1,0,0,1,0,1};
    for(int i = 0; i < s.size(); i++) pass &= exp3[i] == s[i];
    s = msequence(4);
    int exp4[15] = {1,1,1,0,0,0,1,0,0,1,1,0,1,0,1};
    for(int i = 0; i < s.size(); i++) pass &= exp4[i] == s[i];
    s = msequence(5);
    int exp5[31] = {1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,0,0,1,0,0,1,0,1,1,0,0,1};
    for(int i = 0; i < s.size(); i++) pass &= exp5[i] == s[i];
    if(pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void mpilotsequenceTest() {
    std::cout << "pilot m-sequence ... ";
    double tol = 1e-8;
    bool pass = true;
    std::vector<std::complex<double>> pilots = pilotmsequence(3);
    std::complex<double> exp3[7] = {1,1,-1,-1,1,-1,1};
    for(int i = 0; i < pilots.size(); i++) pass &= std::abs(exp3[i] - pilots[i]) < tol;
    if(pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

void innerproductTest() {
    std::cout << "inner product ... ";
    std::vector<std::complex<double>> x = {1,2,3,4};
    std::vector<std::complex<double>> y = {1,std::complex<double>(0,1),1,std::complex<double>(0,1)};
    std::complex<double> expected = std::complex<double>(4,-6);
    bool pass = std::abs(expected - inner_product(x,y)) < 1e-8;
    if(pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL" << std::endl;
}

int main(int argc, char** argv) {

    trapezoidalSquareTest();
    extended_gcd_Test();
    gcdTest();
    fzeroLinearTest();
    fzeroCubicTest();
    msequenceTest();
    mpilotsequenceTest();
    innerproductTest();
    
    return (EXIT_SUCCESS);
}

