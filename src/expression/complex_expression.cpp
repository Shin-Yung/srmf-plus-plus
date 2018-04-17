/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-29 21:51:20
* Last Modified by:   amedhi
* Last Modified time: 2017-05-30 11:49:04
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/

#include "./complex_expression.h"

namespace expr {

void ComplexExpr::add_var(const std::string& name, const double& val) 
{ 
  vars_[name] = mup::Value(val);
  ParserX::DefineVar(name, mup::Variable(&vars_[name]));
}

void ComplexExpr::set_val(const std::string& name, const double& val) 
{ 
  vars_.at(name) = mup::Value(val);
}

void ComplexExpr::set_expr(const std::string& expr)
{
  ParserX::SetExpr(expr);
}

std::complex<double> ComplexExpr::evaluate(void) 
{
  try {
    return ParserX::Eval().GetComplex();
  }
  catch (std::exception& e) 
  { 
    std::string msg = "ComplexExpr::evaluate:\n" + std::string(e.what());
    throw std::runtime_error(msg);
  }
}


} // end namespace expr
