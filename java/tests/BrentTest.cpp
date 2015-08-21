/* 
 * File:   BrentTest.cpp
 * Author: Robby McKilliam
 *
 * Created on 26/01/2013, 3:48:58 PM
 */

#include <stdlib.h>
#include <iostream>
#include "Util.h"

  void QuadraticTest() {
    double tol = 1e-7;
    auto f = [] (double x) { return (x-2.0)*(x-2.0); };
    Brent opt(f, -4.0, 1.0, 5.0);
    std::cout << "Brent QuadraticTest ... ";
    bool pass = fabs(0.0 - opt.fmin()) < tol;
    pass &= fabs(2.0 - opt.xmin()) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL " << opt.fmin() << " " << opt.xmin() << std::endl;
  }
  
    void QuarticTest() {
    double tol = 1e-7;
    auto f = [] (double x) { return x*x*x*x; };
    Brent opt(f, -2.0, 1.0, 2.0);
    std::cout << "Brent QuarticTest ... ";
    bool pass = fabs(0.0 - opt.fmin()) < tol;
    pass &= fabs(0.0 - opt.xmin()) < tol;
    if (pass) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL " << opt.fmin() << " " << opt.xmin() << std::endl;
  }

int main(int argc, char** argv) {
    QuadraticTest();
    QuarticTest();

    return (EXIT_SUCCESS);
}

