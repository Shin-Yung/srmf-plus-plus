/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   amedhi
* Last Modified time: 2017-06-11 16:45:15
*----------------------------------------------------------------------------*/
#include "model.h"

namespace model {


unsigned Hamiltonian::add_siteterm(const std::string& name, const CouplingConstant& cc, 
  const op::quantum_op& op)
{
  // remap site type values in 'cc'
  CouplingConstant cc_remapped = cc;
  cc_remapped.clear_map();
  if (cc.size()==1 && cc.begin()->first==CouplingConstant::global_type) {
    // the 'cc' is implicitly defined for all types
    strMatrix expr_mat = cc.begin()->second;
    for (const auto& m : sitetypes_map_) {
      cc_remapped.insert({m.second, expr_mat});
    }
  }
  else {
    for (auto it=cc.begin(); it!=cc.end(); ++it) {
      unsigned sitetype = it->first;
      auto it2=sitetypes_map_.find(sitetype);
      if (it2!=sitetypes_map_.end()) {
        unsigned mapped_type = it2->second;
        //std::cout << mapped_type << it->second << "\n";
        cc_remapped.insert({mapped_type, it->second});
      }
      else throw std::range_error("Hamiltonian::add_siteterm: non-existent 'site type' specified");
    }
  }
  unsigned num_sitetypes = sitetypes_map_.size();
  this->site_terms_.push_back(HamiltonianTerm(name,op,cc_remapped,num_sitetypes));
  return this->site_terms_.size();
}

unsigned Hamiltonian::add_bondterm(const std::string& name, const CouplingConstant& cc,
  const op::quantum_op& op)
{
  // remap bond type values in 'cc'
  CouplingConstant cc_remapped = cc;
  cc_remapped.clear_map();
  if (cc.size()==1 && cc.begin()->first==CouplingConstant::global_type) {
    // the 'cc' is implicitly defined for all types
    strMatrix expr_mat = cc.begin()->second;
    for (const auto& m : bondtypes_map_) {
      //std::cout << "m= " << m.second << " cc_expr = " << cc_expr << "\n";
      cc_remapped.insert({m.second, expr_mat});
    }
  }
  else {
    for (auto it=cc.begin(); it!=cc.end(); ++it) {
      unsigned bondtype = it->first;
      auto it2=bondtypes_map_.find(bondtype);
      if (it2!=bondtypes_map_.end()) {
        unsigned mapped_type = it2->second;
        cc_remapped.insert({mapped_type, it->second});
      }
      else throw std::range_error("Hamiltonian::add_bondterm: non-existent 'site type' specified");
    }
  }
  unsigned num_bondtypes = bondtypes_map_.size();
  bond_terms_.push_back(HamiltonianTerm(name, op, cc_remapped, num_bondtypes));
  return bond_terms_.size();
}

unsigned Hamiltonian::add_disorder_term(const std::string& name, const op::quantum_op& op)
{
  // site disorder term
  // set dummy 'cc' 
  CouplingConstant cc;
  for (const auto& m : sitetypes_map_) {
    cc.insert({m.second, strMatrix("0.0")});
  }
  unsigned num_sitetypes = sitetypes_map_.size();
  this->disorder_terms_.push_back(HamiltonianTerm(name,op,cc,num_sitetypes));
  // only one disorder term implemented currently
  if (this->disorder_terms_.size()>1)
    throw std::range_error("Hamiltonian::add_disorderterm: only one site disorder term allowed");
  have_disorder_term_ = true;
  return this->disorder_terms_.size();
}

int Hamiltonian::init(const lattice::Lattice& lattice)
{
  // reset
  parms_.clear();
  //operators.clear();
  // maps of site & bond type values (to contigous type values)
  sitetypes_map_ = lattice.sitetypes_map();
  bondtypes_map_ = lattice.bondtypes_map();

  //  sitetypes to sitedim map
  /*type_dim_map_.clear();
  for (unsigned i=0; i<lattice.num_basis_sites(); ++i) {
    unsigned type = lattice.basis_site(i).type();
    unsigned dim = lattice.basis_site(i).num_orbitals();
    unsigned mapped_type = sitetypes_map_.at(type);
    type_dim_map_.insert({mapped_type, dim});
  }*/

  // maps of a given bond type to the types of its target
  /*bond_sites_map_.clear();
  for (unsigned i=0; i<lattice.num_basis_bonds(); ++i) {
    lattice::Bond b = lattice.basis_bond(i);
    lattice::Site src = lattice.basis_site(b.src_id());
    lattice::Site tgt = lattice.basis_site(b.tgt_id());
    bond_sites_map_.insert({b.type(), std::make_pair(src.type(), tgt.type())});
    //std::cout << "bond_site_map = "<<b.type()<<" "<<src.type()<<" "<<tgt.type()<<"\n";
  }*/
  bond_terms_.clear();
  site_terms_.clear();
  disorder_terms_.clear();
  bt_begin_ = bond_terms_.cbegin();
  bt_end_ = bond_terms_.cend();
  st_begin_ = site_terms_.cbegin();
  st_end_ = site_terms_.cend();
  dterm_begin_ = disorder_terms_.cbegin();
  dterm_end_ = disorder_terms_.cend();
  return 0;
}

int Hamiltonian::finalize(const lattice::Lattice& L)
{
  // check if 'sitebasis' for all 'site types' are defined
  /*for (const auto& elem : sitetypes_map_) {
    unsigned site_type = elem.second;
    if (basis_.find(site_type) == basis_.end()) 
      throw std::range_error("modellibrary:: 'sitebasis' not defined for all 'site type'-s");
  }*/
  // finalize the site terms
  for (auto it=site_terms_.begin(); it!=site_terms_.end(); ++it) {
    it->eval_coupling_constant(constants_, parms_); 
  }
  have_siteterm_ = (site_terms_.size()>0);
  st_begin_ = site_terms_.cbegin();
  st_end_ = site_terms_.cend();
  // finalize the bond terms
  for (auto it=bond_terms_.begin(); it!=bond_terms_.end(); ++it) {
    it->eval_coupling_constant(constants_, parms_); 
  }
  have_bondterm_ = (bond_terms_.size()>0);
  bt_begin_ = bond_terms_.cbegin();
  bt_end_ = bond_terms_.cend();

  // disorder terms
  dterm_begin_ = disorder_terms_.cbegin();
  dterm_end_ = disorder_terms_.cend();

  // info string
  set_info_string(L);
  return 0;
}

/*void Hamiltonian::change_parameter_value(const std::string& pname, const double& pval) 
{
  auto it = parms_.find(pname);
  if (it != parms_.end()) it->second = pval;
}*/

void Hamiltonian::update_parameter(const std::string& pname, const double& val)
{
  parms_.at(pname) = val; 
}

void Hamiltonian::update_parameters(const input::Parameters& inputs)
{
  // update the parameter values
  int info;
  for (auto& p : parms_) p.second = inputs.set_value(p.first, p.second, info);
}

void Hamiltonian::update_terms(void)
{
  // update the model term couping constants
  for (auto it=site_terms_.begin(); it!=site_terms_.end(); ++it) {
    it->eval_coupling_constant(constants_, parms_); 
  }
  for (auto it=bond_terms_.begin(); it!=bond_terms_.end(); ++it) {
    it->eval_coupling_constant(constants_, parms_); 
  }
}

double Hamiltonian::get_parameter_value(const std::string& pname) const
{
  auto it = parms_.find(pname);
  if (it != parms_.end()) return it->second;
  else throw std::logic_error("Hamiltonian::get_parameter_value: parameter does not exist");
}

void Hamiltonian::get_term_names(std::vector<std::string>& term_names) const
{
  term_names.clear();
  for (auto it=bond_terms_.cbegin(); it!=bond_terms_.cend(); ++it) 
    term_names.push_back(it->name());
  for (auto it=site_terms_.cbegin(); it!=site_terms_.cend(); ++it) 
    term_names.push_back(it->name());
  for (auto it=disorder_terms_.cbegin(); it!=disorder_terms_.cend(); ++it) 
    term_names.push_back(it->name());
}

void Hamiltonian::set_info_string(const lattice::Lattice& L) 
{
  std::ostringstream info_ostream;
  info_ostream.clear();
  info_ostream << "# Lattice: " << L.name() << " (";
  info_ostream << "Size="<<L.size1()<<"x"<<L.size2()<<"x"<< L.size3()<<", ";
  info_ostream << "Sites/unitcell="<<L.num_basis_sites()<<", ";
  info_ostream << "Boundary="<<static_cast<int>(L.bc1_periodicity()) << "-"; 
  info_ostream << static_cast<int>(L.bc2_periodicity()) << "-";
  info_ostream << static_cast<int>(L.bc3_periodicity()) << ")\n";
  info_ostream << "# No of sites = " << L.num_sites() << "\n";
  info_ostream << "# Model: " << model_name << "\n";
  info_ostream.precision(6);
  info_ostream.setf(std::ios_base::fixed);
  for (const auto& p : parms_) 
    info_ostream << "# " << p.first << " = " << p.second << "\n";
  info_str_ = info_ostream.str();
  // signature string
  std::ostringstream signature_ostream;
  signature_ostream << "L" << static_cast<int>(L.id()) << "_"; 
  signature_ostream << L.size1() << "x" << L.size2() << "x" << L.size3();
  signature_ostream << "_" << model_name;
  signature_ostream.precision(3);
  signature_ostream.setf(std::ios_base::fixed);
  for (const auto& p : parms_) signature_ostream << "_" << p.first << p.second;
  signature_str_ = signature_ostream.str();
}





} // end namespace model
