/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-16 16:47:10
*----------------------------------------------------------------------------*/
#include "qn.h"

namespace model {

unsigned QuantumNumber::count_ = 0;

QuantumNumber::QuantumNumber(const std::string& name, 
  const value_type& min, const value_type& max, const value_type& step, const bool& fermionic)
  : name_{name}, min_{min}, max_{max}, step_{step}, fermionic_{fermionic}
{
  if (step<0 || min_>max_)  
    throw std::invalid_argument("QuantumNumber::QuantumNumber");
  if (step_==0 && min_ != max_) 
    throw std::invalid_argument("QuantumNumber::QuantumNumber");
  if (max_>min_ && (max_-min_)%step_ !=0)
    throw std::invalid_argument("QuantumNumber::QuantumNumber");

  if (min_ == max_) num_states_ = 1;
  else num_states_ = (max_-min_)/step_ + 1;
  id_ = count_++;
}


} // end namespace model
