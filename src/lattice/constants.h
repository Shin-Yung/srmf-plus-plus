/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-01-26 10:23:16
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-15 11:44:48
*----------------------------------------------------------------------------*/
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <complex>

const double dp_tol = 1.0E-12;
//constexpr double pi() { return 4.0 * std::atan(1.0); }
//constexpr double two_pi() { return 8.0 * std::atan(1.0); }
//constexpr double half_pi() { return 2.0 * std::atan(1.0); }
constexpr double pi(void) { return 3.1415926535897932384626433832795028841971693993751058209; }
constexpr double two_pi(void) { return 2.0 * pi(); }
constexpr double half_pi(void) { return 0.5 * pi(); }
//constexpr std::complex<double> ii(void) { return std::complex<double>(0.0, 1.0); }


#endif
