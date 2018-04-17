/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-16 22:12:57
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-27 23:49:33
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
//#define EIGEN_USE_MKL_ALL
#include <Eigen/Core>

#ifndef MATRIX_H
#define MATRIX_H

#define REAL_WAVEFUNCTION

using Vector3i = Eigen::Vector3i;
using Vector3d = Eigen::Vector3d;
using Diagonal3d = Eigen::DiagonalMatrix<double,3>;
using RealVector = Eigen::VectorXd;
using RealMatrix = Eigen::MatrixXd;
using ComplexVector = Eigen::VectorXcd;
using ComplexMatrix = Eigen::MatrixXcd;

#ifdef REAL_WAVEFUNCTION
  using amplitude_t = double;
  #define ampl_part(psi) std::real(psi) 
  using Matrix = Eigen::MatrixXd;
  using Vector = Eigen::VectorXd;
  using ColVector = Eigen::VectorXd;
  using RowVector = Eigen::RowVectorXd;
#else
  using amplitude_t = std::complex<double>;
  #define ampl_part(psi) psi 
  using Matrix = Eigen::MatrixXcd;
  using Vector = Eigen::VectorXcd;
  using ColVector = Eigen::VectorXcd;
  using RowVector = Eigen::RowVectorXcd;
#endif

#endif
