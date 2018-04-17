/*---------------------------------------------------------------------------
* Parameters: Handles the parameter values for individual tasks.
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-27 00:31:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-09 00:33:04
*----------------------------------------------------------------------------*/
// File: taskparams.cc

#include <cmath>
#include "inputparams.h"

namespace input {

int Parameters::set_value(const std::string& pname, const int& defval) const
{
  it = this->find(pname);
  if (it != this->end()) {
    switch (it->second.type) {
      case ptype::num: return static_cast<int>(round(it->second.num_val));
      case ptype::boo: return it->second.bool_val;
      default: warn_type_mismatch(pname, "int or bool"); return defval;
    }
  }
  warn_not_found(pname); return defval;
}

int Parameters::set_value(const std::string& pname, const int& defval, int& info) const
{
  info = 0;
  it = this->find(pname);
  if (it != this->end()) {
    switch (it->second.type) {
      case ptype::num: return static_cast<int>(round(it->second.num_val));
      case ptype::boo: return it->second.bool_val;
      default: info = 1; return defval;
    }
  }
  info = 1; return defval;
}

int Parameters::set_value(const std::string& pname, const int& defval, int& info, bool& is_const) const
{
  info = 0;
  it = this->find(pname);
  if (it != this->end()) {
    is_const = it->second.is_const;
    switch (it->second.type) {
      case ptype::num: return static_cast<int>(round(it->second.num_val));
      case ptype::boo: return it->second.bool_val;
      default: info = 1; return defval;
    }
  }
  is_const = false; info = 1; return defval;
}

double Parameters::set_value(const std::string& pname, const double& defval) const
{
  it = this->find(pname);
  if (it != this->end()) {
    if (it->second.type == ptype::num) return it->second.num_val;
    warn_type_mismatch(pname, "double"); return defval;
  }
  warn_not_found(pname); return defval; 
}

double Parameters::set_value(const std::string& pname, const double& defval, int& info) const
{
  it = this->find(pname);
  if (it != this->end() && it->second.type == ptype::num) {
    info = 0; return it->second.num_val;
  }
  info = 1; return defval; 
}

double Parameters::set_value(const std::string& pname, const double& defval, int& info, 
  bool& is_const) const
{
  it = this->find(pname);
  if (it != this->end() && it->second.type == ptype::num) {
    is_const = it->second.is_const; 
    info = 0; return it->second.num_val;
  }
  is_const = false; info = 1; return defval; 
}

std::string Parameters::set_value(const std::string& pname, const std::string& defval) const
{
  it = this->find(pname);
  if (it != this->end()) {
    if (it->second.type == ptype::str) return it->second.str_val;
    warn_type_mismatch(pname, "std::string"); return defval;
  }
  warn_not_found(pname); return defval; return defval;
}

std::string Parameters::set_value(const std::string& pname, const std::string& defval, int& info) const
{
  it = this->find(pname);
  if (it != this->end()) {
    if (it->second.type == ptype::str) return it->second.str_val;
    warn_type_mismatch(pname, "std::string"); return defval;
  }
  info = 1; return defval; 
}

bool Parameters::is_constant(const std::string& pname) const
{
  it = this->find(pname);
  if (it != this->end()) return it->second.is_const;
  warn_not_found(pname); return false;
}

void Parameters::warn_not_found(const std::string& pname) const 
{
  std::cout << ">> alert: Parameters::set_value: parameter '";
  std::cout << pname << "' not found, default value chosen\n";
}

void Parameters::warn_type_mismatch(const std::string& pname, const std::string& type) const 
{
  std::cout << ">> alert: Parameters::set_value: parameter '";
  std::cout << pname << "' of type '"  << type << "' not found, default value chosen\n";
}

void Parameters::show(const unsigned& set_id) const
{
  std::cout << "Parameter set = " << set_id+1 << ":" << std::endl; 
  std::cout << "---------------------" << std::endl; 
  std::map<std::string, pval>::const_iterator it;
  for (it=this->begin(); it!=this->end(); ++it) {
    std::cout << " " << it->first << " = "; 
    switch (it->second.type) {
      case ptype::boo:
        std::cout << it->second.bool_val; break;
      case ptype::num:
        std::cout << it->second.num_val; break;
      case ptype::str:
        std::cout << "\"" << it->second.str_val << "\""; break;
      case ptype::nan:
        throw std::logic_error("Undefined parameter type detected"); 
        break;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
}


} // end namespace input
