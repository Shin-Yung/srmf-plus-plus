/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef BASIS_H
#define BASIS_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "qn.h"
#include "quantum_operator.h"

namespace model {

class SiteState : public std::vector<QuantumNumber::value_type> {
public:
  using value_type = QuantumNumber::value_type;
  SiteState() {}
  SiteState(const std::vector<value_type>& x) : std::vector<value_type>(x)  {}
  ~SiteState() {}
};

class SiteBasis : public std::vector<QuantumNumber>
{
public:
  using QN = QuantumNumber;
  using value_type = QN::value_type;
  SiteBasis() : num_states_{0}, valid_{false} {}
  SiteBasis(const std::string& name, const value_type& min, const value_type& max, 
    const value_type& step, const bool& fermionic=false);
  ~SiteBasis() {}
  void create(const std::string& name, const value_type& min, const value_type& max, 
    const value_type& step, const bool& fermionic=false);
  unsigned add_qn(const std::string& name, const value_type& min, const value_type& max, 
    const value_type& step, const bool& fermionic=false);
  unsigned add_qn(const QN& qn);
  unsigned add_operator(const std::string& name, const std::string& matrixelem_expr, 
    const std::string& qn, const QN::value_type& change=0);
  void clear(void);
  double apply(const unsigned& op_id, const unsigned& state_idx) const;
  double apply(const std::string& op, const unsigned& state_idx) const;
  const std::string& qn_name(const unsigned& i) const { return at(i).name(); }
  const std::string& op_name(const unsigned& i) const { return operators_.at(i).name(); }
  const unsigned& dimension(void) const { return num_states_; }
  unsigned num_operator(void) const { return operators_.size(); }
  const QuantumOperator& quantum_operator(const unsigned& i) const 
  { return operators_[i]; }
  const bool& valid(void) const { return valid_; }
  void finalize(void);
private:
  unsigned num_states_{0};
  bool valid_{false};
  using op_result = std::pair<double, QN::value_type>; // (matrix_element, new_state) pair
  std::vector<QuantumOperator> operators_;
  std::map<unsigned, std::vector<op_result> > operator_action_;
  std::vector<SiteState> basis_states_;
  //operator_map operators_;
};

class BasisDescriptor : public std::map<unsigned,SiteBasis>
{
public:
  //using QuantumNumber = QuantumNitumber;
  BasisDescriptor() : name_{""} {}
  ~BasisDescriptor() {}
  bool add_sitebasis(const unsigned& type, SiteBasis& sitebasis)
  {
    sitebasis.finalize();
    auto res = insert(std::make_pair(type, sitebasis));
    //return this->size();
    return res.second;
  }
  const unsigned& dimension(const unsigned& site_type) const 
  {
    auto it = find(site_type);
    if (it != end())
      return it->second.dimension();
    else throw std::runtime_error("BasisDescriptor::dimension: given 'site_type' not found");
  }
  const SiteBasis& site_basis(const unsigned& type) const 
  { return this->at(type); }
private:
  std::string name_;
};



} // end namespace model

#endif
