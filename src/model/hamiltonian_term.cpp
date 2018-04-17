/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   amedhi
* Last Modified time: 2017-05-30 11:59:41
*----------------------------------------------------------------------------*/
#include "hamiltonian_term.h"
//#include "../expression/expression.h"
#include "../expression/complex_expression.h"

namespace model {

/*--------------------------CouplingConstant--------------------*/
const int CouplingConstant::global_type = -1;

CouplingConstant::CouplingConstant(const std::string& expr)
{
  super_type::clear();
  // expr is applicable for all site & bond types
  super_type::insert({global_type, strMatrix(expr)});
  num_types_ = -1;
  valid_ = true;
}

CouplingConstant::CouplingConstant(const strMatrix::row_t& expr_vec)
{
  super_type::clear();
  // expr is applicable for all site & bond types
  super_type::insert({global_type, strMatrix(expr_vec)});
  num_types_ = -1;
  valid_ = true;
}

CouplingConstant& CouplingConstant::operator=(const std::string& expr)
{
  super_type::clear();
  // expr is applicable for all site & bond types
  super_type::insert({global_type, strMatrix(expr)});
  num_types_ = -1; 
  valid_ = true;
  return *this;
}

CouplingConstant& CouplingConstant::operator=(const strMatrix::row_t& expr_vec)
{
  super_type::clear();
  // expr is applicable for all site & bond types
  super_type::insert({global_type, strMatrix(expr_vec)});
  num_types_ = -1; 
  valid_ = true;
  return *this;
}

CouplingConstant& CouplingConstant::operator=(const strMatrix& expr_mat)
{
  super_type::clear();
  // expr is applicable for all site & bond types
  super_type::insert({global_type, expr_mat});
  num_types_ = -1; 
  valid_ = true;
  return *this;
}

void CouplingConstant::create(const unsigned& num_types) 
{
  super_type::clear();
  num_types_ = num_types;
  valid_ = false;
}

void CouplingConstant::add_type(const unsigned& type, const std::string& expr)
{
  super_type::insert({type, strMatrix(expr)});
  valid_ = (num_types_==static_cast<int>(size()));
}

void CouplingConstant::add_type(const unsigned& type, 
  const std::vector<std::string>& expr_vec)
{
  super_type::insert({type, {expr_vec}});
  valid_ = (num_types_==static_cast<int>(size()));
}

void CouplingConstant::add_type(const unsigned& type, 
  const strMatrix& expr_mat)
{
  super_type::insert({type, expr_mat});
  valid_ = (num_types_==static_cast<int>(size()));
}


//-----------------------HamiltonianTerm-------------------------
void HamiltonianTerm::construct(const std::string& name, const op::quantum_op& op, 
    const CouplingConstant& cc, const unsigned& size)
{
  if (!cc.valid()) throw std::invalid_argument("HamiltonianTerm:: Invalid CouplingConstant");
  name_ = name;
  op_ = op;
  cc_ = cc;
  max_operand_types_ = size;

  cc_values_.resize(max_operand_types_);
  is_defined_.resize(max_operand_types_);
  for (unsigned i=0; i<max_operand_types_; ++i) {
    is_defined_[i] = false;
  }

  // if the 'cc' is implicitly defined for all types 
  if (cc_.size()==1 && cc_.begin()->first==CouplingConstant::global_type) {
    int rows = cc_.begin()->second.rows();
    int cols = cc_.begin()->second.cols();
    for (unsigned i=0; i<max_operand_types_; ++i) {
      cc_values_[i] = ComplexMatrix::Zero(rows, cols);
      is_defined_[i] = true;
    }
  } 
  else {
    // operator is defined only for those types for which 'cc' is set explicitly
    for (const auto& p : cc_) {
      cc_values_[p.first] = ComplexMatrix::Zero(p.second.rows(), p.second.cols());
      is_defined_[p.first] = true;
    }
  }
}

void HamiltonianTerm::eval_coupling_constant(const ModelParams& cvals, const ModelParams& pvals)
{
  expr::ComplexExpr expr;
  for (const auto& c : cvals) expr.add_var(c.first, c.second);
  for (const auto& p : pvals) expr.add_var(p.first, p.second);
  try { 
    // if the 'cc' is implicitly defined for all types 
    if (cc_.size()==1 && cc_.begin()->first==CouplingConstant::global_type) {
      int rows = cc_.begin()->second.rows();
      int cols = cc_.begin()->second.cols();
      ccval_t mat(rows, cols);
      for (int i=0; i<rows; ++i) {
        for (int j=0; j<cols; ++j) {
          expr.set_expr(cc_.begin()->second(i,j));
          auto val = expr.evaluate();
          mat(i,j) = val;
        }
      }
      for (auto& v : cc_values_) v = mat;
    }
    else {
      for (const auto& p : cc_) {
        int rows = p.second.rows();
        int cols = p.second.cols();
        ccval_t mat(rows, cols);
        for (int i=0; i<rows; ++i) {
          for (int j=0; j<cols; ++j) {
            expr.set_expr(p.second(i,j));
            auto val = expr.evaluate();
            mat(i,j) = val;
          }
        }
        cc_values_[p.first] = mat;
      }
    }
  }
  catch (std::exception& e) 
  { 
    std::string msg = "BondOperatorTerm::evaluate_coupling_constant:\n" + std::string(e.what());
    throw std::runtime_error(msg);
  }
}



} // end namespace model
