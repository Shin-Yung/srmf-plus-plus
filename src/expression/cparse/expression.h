/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef EXPRESSION_H
#define EXPRESSION_H
#include <iostream>
#include <string>
#include "shunting_yard.h"

namespace expr {

class Expression : public calculator
{
public:
  using variables = TokenMap; 
  Expression() {}
  Expression(const std::string& expr_str) : calculator(expr_str.c_str()) {}
  ~Expression() {}
  void set_variable(const std::string& x, const double& val) 
  { vars_[x] = val; }
  double evaluate(const std::string& expr_str) 
  { 
    try {
      return this->calculator::calculate(expr_str.c_str(), &vars_).asDouble(); 
    }
    catch (std::exception& e) {
      std::string msg = "Expression Error: " + std::string(e.what());
      throw std::runtime_error(msg);
    }
  }
  
  double evaluate(const std::string& expr_str, TokenMap& vars) 
  { 
    try {
      return this->calculator::calculate(expr_str.c_str(), &vars).asDouble(); 
    }
    catch (std::exception& e) {
      std::string msg = "Expression Error: " + std::string(e.what());
      throw std::runtime_error(msg);
    }
  }

  double evaluate(TokenMap& var) 
  { 
    try {
      return this->calculator::eval(&var).asDouble(); 
    }
    catch (...) {
      std::cout << "Expression Error: "; throw;
    }
  }

  int evaluate_asInt(TokenMap& var) 
  { 
    try {
      return this->calculator::eval(&var).asInt(); 
    }
    catch (...) {
      std::cout << "Expression Error: "; throw;
    }
  }
private:
  TokenMap vars_;
};

} // end namespace expr

#endif
