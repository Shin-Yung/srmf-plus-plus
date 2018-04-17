/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef QN_H
#define QN_H

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

namespace model {

class QuantumNumber 
{
public:
  using value_type = short;
  QuantumNumber(const std::string& name, const value_type& min=0, 
    const value_type& max=0, const value_type& step=1, const bool& fermionic=false);
  ~QuantumNumber() {};
  const value_type& min(void) const { return min_; }
  const value_type& max(void) const { return max_; }
  const value_type& step(void) const { return step_; }
  const unsigned& id(void) const { return id_; }
  const std::string& name(void) const { return name_; }
  const std::size_t& num_states(void) const { return num_states_; }
  const bool& fermionic(void) const { return fermionic_; }
private:
  static unsigned count_;
  unsigned id_;
  std::string name_;
  value_type min_;
  value_type max_;
  value_type step_;
  bool fermionic_;
  //bool valid_;
  std::size_t num_states_;
};


} // end namespace model

#endif
