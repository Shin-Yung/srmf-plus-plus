/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:46
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-14 11:01:07
*----------------------------------------------------------------------------*/
#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <sstream>
#include <string>
#include <complex>
#include <vector>
#include <map>
#include <stdexcept>
#include "strmatrix.h"
#include "../lattice/lattice.h"
#include "./hamiltonian_term.h"

namespace model {

enum class model_id {
  UNDEFINED, HUBBARD, TJ, DISORDERED_TJ, TBI_HUBBARD
};

class Hamiltonian 
{
public:
  using siteterm_iterator = std::vector<HamiltonianTerm>::const_iterator; 
  using bondterm_iterator = std::vector<HamiltonianTerm>::const_iterator; 
  Hamiltonian() {}
  Hamiltonian(const input::Parameters& inputs, const lattice::Lattice& lattice)
  { construct(inputs, lattice); }
  ~Hamiltonian() {}
  int construct(const input::Parameters& inputs, const lattice::Lattice& lattice);
  virtual int init(const lattice::Lattice& lattice);
  int define_model(const input::Parameters& inputs, const lattice::Lattice& lattice);
  int finalize(const lattice::Lattice& lattice);

  //unsigned add_sitebasis(SiteBasis& sitebasis);
  //unsigned add_sitebasis(const unsigned& type, SiteBasis& sitebasis);
  unsigned add_parameter(const std::string& pname, const double& defval, 
    const input::Parameters& inputs)
    { parms_[pname] = inputs.set_value(pname, defval); return parms_.size(); }
  unsigned add_parameter(const std::string& pname, const double& defval, 
    const input::Parameters& inputs, int& info)
    { parms_[pname] = inputs.set_value(pname, defval, info); return parms_.size(); }
  unsigned add_parameter(const std::string& pname, const double& val) 
    { parms_[pname] = val; return parms_.size(); }
  void update_parameters(const input::Parameters& inputs);
  void update_parameter(const std::string& pname, const double& val); 
  virtual void update_terms(void);
  //void change_parameter_value(const std::string& pname, const double& pval);
  double get_parameter_value(const std::string& pname) const;
  unsigned add_constant(const std::string& cname, const double& val) 
    { constants_.insert({cname, val}); return constants_.size(); }
  unsigned add_siteterm(const std::string& name, const CouplingConstant& cc, const op::quantum_op& op);
  unsigned add_bondterm(const std::string& name, const CouplingConstant& cc, const op::quantum_op& op);
  unsigned add_disorder_term(const std::string& name, const op::quantum_op& op);
  void set_no_dbloccupancy(void) { double_occupancy_=false; }

  //const BasisDescriptor& basis(void) const { return basis_; }
  //const SiteBasis& site_basis(const unsigned& site_type) const { return basis_.at(site_type); }
  //unsigned sitebasis_dimension(const unsigned& site_type) const
  //{ return basis_.dimension(site_type); }
  const ModelParams& parameters(void) const { return parms_; }
  const ModelParams& constants(void) const { return constants_; }
  const bool& double_occupancy(void) const { return double_occupancy_; }
  const bool& have_disorder_term(void) const { return have_disorder_term_; }
  const bool& have_siteterm(void) const { return have_siteterm_; }
  const bool& have_bondterm(void) const { return have_bondterm_; }
  const siteterm_iterator& siteterms_begin(void) const { return st_begin_; }
  const siteterm_iterator& siteterms_end(void) const { return st_end_; }
  const bondterm_iterator& bondterms_begin(void) const { return bt_begin_; }
  const bondterm_iterator& bondterms_end(void) const { return bt_end_; }
  const siteterm_iterator& disorder_term_begin(void) const { return dterm_begin_; }
  const siteterm_iterator& disorder_term_end(void) const { return dterm_end_; }
  std::pair<siteterm_iterator, siteterm_iterator> site_terms(void) const 
    { return std::make_pair(site_terms_.cbegin(), site_terms_.cend()); }
  std::pair<bondterm_iterator, bondterm_iterator> bond_terms(void) const 
    { return std::make_pair(bond_terms_.cbegin(), bond_terms_.cend()); }
  std::pair<siteterm_iterator, siteterm_iterator> disorder_terms(void) const 
    { return std::make_pair(disorder_terms_.cbegin(), disorder_terms_.cend()); }
  unsigned num_siteterms(void) const { return site_terms_.size(); }
  unsigned num_bondterms(void) const { return bond_terms_.size(); }
  unsigned num_disorder_terms(void) const { return disorder_terms_.size(); }
  unsigned num_terms(void) const 
    { return site_terms_.size()+bond_terms_.size()+disorder_terms_.size(); }
  void get_term_names(std::vector<std::string>& term_names) const;
  const std::string& signature_str(void) const { return signature_str_; }
  std::ostream& print_info(std::ostream& os) const { return os << info_str_; }
  //const BondTerm::BondSiteMap& bond_sites_map(void) const { return bond_sites_map_; }

  //const SiteTerm& siteterm(const unsigned& i) const { return siteterms_[i]; };
  //const BondTerm& bondterm(const unsigned& i) const { return bondterms_[i]; };
private:
  model_id mid {model_id::UNDEFINED};
  std::string model_name;
  //BasisDescriptor basis_;
  std::map<unsigned, unsigned> sitetypes_map_;
  std::map<unsigned, unsigned> bondtypes_map_;
  std::map<unsigned, unsigned> type_dim_map_;
  //BondTerm::BondSiteMap bond_sites_map_;  
  std::vector<HamiltonianTerm> bond_terms_;
  std::vector<HamiltonianTerm> site_terms_;

  bool double_occupancy_{true};
  bool have_siteterm_{false};
  bool have_bondterm_{false};
  siteterm_iterator st_begin_;
  siteterm_iterator st_end_;
  bondterm_iterator bt_begin_;
  bondterm_iterator bt_end_;
  //std::vector<SiteTermDescriptor> siteterms_;
  //std::vector<BondTermDescriptor> bondterms_;
  ModelParams parms_;
  ModelParams constants_;

  // disorder term
  bool have_disorder_term_{false};
  //SiteDisorder site_disorder_;
  std::vector<HamiltonianTerm> disorder_terms_;
  siteterm_iterator dterm_begin_;
  siteterm_iterator dterm_end_;

  //std::ostringstream info_str_;
  //std::ostringstream signature_str_;
  std::string info_str_;
  std::string signature_str_;
  void set_info_string(const lattice::Lattice& lattice); 
};


} // end namespace model

#endif