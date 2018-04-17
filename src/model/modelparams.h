/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:46
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-16 16:47:10
*----------------------------------------------------------------------------*/
#ifndef MODELPARAMS_H
#define MODELPARAMS_H

#include <iostream>
#include <string>
#include <unordered_map>
//#include <complex>

namespace model {

class ModelParams : public std::unordered_map<std::string, double>
{
public:
  using super_type = std::unordered_map<std::string, double>;
  using iterator = super_type::iterator;
  using const_iterator = super_type::const_iterator;
  ModelParams() {}
  ~ModelParams() {}
private:
};

} // end namespace model

#endif
