/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-29 21:48:38
* Last Modified by:   amedhi
* Last Modified time: 2017-05-30 11:49:11
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef COMPLEX_EXPRESSION_H
#define COMPLEX_EXPRESSION_H
#include <iostream>
#include <string>
#include <complex>
#include "./muparserx/parser/mpParser.h"

namespace expr {

class ComplexExpr : public mup::ParserX
{
public:
  using variables = std::map<std::string,mup::Value>; 
  ComplexExpr() : ParserX(mup::pckALL_COMPLEX) {}
  ~ComplexExpr() {}
  void add_var(const std::string& name, const double& val); 
  void set_val(const std::string& name, const double& val); 
  void set_expr(const std::string& expr);
  std::complex<double> evaluate(void); 
private:
  variables vars_;
};


} // end namespace expr

#endif