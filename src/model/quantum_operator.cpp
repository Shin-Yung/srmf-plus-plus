/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-16 16:47:10
*----------------------------------------------------------------------------*/
#include "quantum_operator.h"
#include "../expression/expression.h"

namespace model {

unsigned QuantumOperator::count_ = 0;

QuantumOperator::QuantumOperator(const std::string& name, const std::string& operand_name)
  : name_(name), matrixelem_(operand_name) 
{
  QuantumNumber::value_type change = 0;
  id_ = count_++;
  this->push_back(std::make_pair(operand_name,change));
}

QuantumOperator::QuantumOperator(const std::string& name, const std::string& elem, 
  const std::string& operand_name, const QN::value_type& change) 
    : name_(name), matrixelem_(elem) 
{ 
  id_ = count_++;
  this->push_back(std::make_pair(operand_name,change));
}

double QuantumOperator::matrix_element(const std::string& qn_name, const QN::value_type& qn_state)
{
  for (unsigned i=0; i<size(); ++i) {
    if (qn_name == operator[](i).first) {
      // operator acts on this quantum number
      expr::Expression matrixelem_expr; //(matrixelem_);
      //expr::Expression::variables vars;
      matrixelem_expr.set_variable(qn_name, qn_state);
      //vars[qn_name] = static_cast<double>(qn_state);
      try {
        //double m = matrixelem_expr.evaluate(vars);
        double m = matrixelem_expr.evaluate(matrixelem_);
        return m; 
      }
      catch (std::exception& e) {
        std::string msg = "QuantumOperator::matrix_element:\n" + std::string(e.what());
        throw std::runtime_error(msg);
      }
    }
  } 
  return 0.0;
}



} // end namespace model
