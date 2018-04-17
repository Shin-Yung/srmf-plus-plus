/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-16 16:47:10
*----------------------------------------------------------------------------*/
#include <tuple>
#include <algorithm>
#include <stdexcept>
#include "sitebasis.h"

namespace model {

SiteBasis::SiteBasis(const std::string& name, const value_type& min, const value_type& max, 
  const value_type& step, const bool& fermionic) 
    : std::vector<QuantumNumber>({QuantumNumber(name,min,max,step,fermionic)})
{
}

void SiteBasis::clear(void)
{
  std::vector<QuantumNumber>::clear();
  num_states_ = 0;
  valid_ = false;
  operators_.clear();
  operator_action_.clear();
  basis_states_.clear();
}

void SiteBasis::create(const std::string& name, const value_type& min, 
  const value_type& max, const value_type& step, const bool& fermionic)
{
  SiteBasis::clear();
  push_back(QuantumNumber(name,min,max,step,fermionic));
}
  
unsigned SiteBasis::add_qn(const std::string& name, const value_type& min, 
  const value_type& max, const value_type& step, const bool& fermionic)
{
  // check if the qunatum number already exist in the basis
  for (const auto& qn : *this) {
    if (qn.name()==name) {
      throw std::range_error("SiteBasis::add_qn: quantum number '"+name+"' already exist"); 
    }
  }
  push_back(QuantumNumber(name,min,max,step,fermionic));
  return back().id();
}
  
unsigned SiteBasis::add_qn(const QuantumNumber& qn)
{
  // check if the qunatum number already exist in the basis
  for (const auto& elem : *this) {
    if (elem.name()==qn.name()) {
      throw std::range_error("SiteBasis::add_qn: quantum number '"+qn.name()+"' already exist"); 
    }
  }
  push_back(qn);
  return back().id();
}

unsigned SiteBasis::add_operator(const std::string& name, const std::string& matrixelem_expr, 
  const std::string& qn, const QN::value_type& change)
{
  operators_.push_back(QuantumOperator(name, matrixelem_expr, qn, change));
  return operators_.back().id();
}

void SiteBasis::finalize(void)
{
  // store all the site basis states
  basis_states_.clear();
  num_states_ = 1;
  for (unsigned i=0; i<size(); ++i) 
    num_states_ *= operator[](i).num_states();
  
  SiteState state;
  for (unsigned i=0; i<size(); ++i) {
    state.push_back(operator[](i).min());
    //std::cout << state[i] << "\n";
  }
  basis_states_.push_back(state);

  unsigned n = 1;
  while (n < num_states_) {
    for (unsigned i=0; i<size(); ++i) {
      if (state[i] < operator[](i).max()) {
        state[i] += operator[](i).step();
        break;
      }
      if (state[i] == operator[](i).max()) state[i] = operator[](i).min();
    }
    basis_states_.push_back(state);
    n++;
  }
  // check the states
  /*for (unsigned i=0; i<basis_states_.size(); ++i) {
    std::cout << "state " << i << ": ";
    for (unsigned j=0; j<size(); ++j) {
      std::cout << basis_states_[i][j] << "  ";
    }
    std::cout << "\n";
  }*/

  //int n = size();
  value_type change;
  std::string qn_name;
  std::vector<op_result> op_results(num_states_);
  for (unsigned i=0; i<operators_.size(); ++i) {
    // matrix element expression string
    // for the time being, assuming one operand (id=0) only
    std::tie(qn_name,change) = operators_[i].action(0);
    // if the 'qn' exist in the basis
    unsigned qn;
    for (qn=0; qn<size(); ++qn) {
      if (qn_name == this->operator[](qn).name()) {
        for (unsigned s=0; s<num_states_; ++s) {
          QN::value_type qn_value = basis_states_[s][qn];
          double matrixelem = operators_[i].matrix_element(qn_name, qn_value);
          // right now assuming only diagonal operators
          op_results[s] = std::make_pair(matrixelem, s);
          //std::cout << qn_name << " " << s << " " << qn_value << "\n";
        }
        // insert the action results
        operator_action_.insert(std::make_pair(operators_[i].id(),op_results));
        break;
      }
    }
    // if 'qn' not found
    if (qn >= size()) 
      throw std::range_error("SiteBasis::finalize: quantum number for operator '"
        +operators_[i].name()+"' not found");
  }
  valid_ = true;
}

double SiteBasis::apply(const unsigned& op_id, const unsigned& state_idx) const
{ 
  return operator_action_.at(op_id).at(state_idx).first;
}

double SiteBasis::apply(const std::string& op_name, const unsigned& state_idx) const
{ 
  for (const auto& op : operators_) {
    if (op.name()==op_name) {
      unsigned op_id = op.id(); 
      //std::cout << operator_action_.at(op_id).size() << "\n";
      //return 1;
      return operator_action_.at(op_id).at(state_idx).first;
    }
  }
  throw std::invalid_argument("SiteBasis::apply Operator '"+op_name+"' not defined for this basis");
}




} // end namespace basis
