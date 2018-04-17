/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef HAMILTONIAN_TERM_H
#define HAMILTONIAN_TERM_H

#include <string>
#include <vector>
#include <complex>
#include <array>
#include <map>
//#include <unordered_map>
#include <stdexcept>
#include <Eigen/Core>
#include "modelparams.h"
#include "strmatrix.h"
#include "../lattice/lattice.h"
#include "../diag/matrix.h"
#include "./quantum_op.h"

namespace model {

enum class qn_op {
  ni_up, ni_dn, ni, cdagicj_up, cdagicj_dn, sisj, niup_nidn, cdagiup_cdagidn, 
  cdagiup_cdagjdn, null
};

class CouplingConstant : public std::unordered_map<int, strMatrix>
{
public:
  using super_type = std::unordered_map<int,strMatrix>;
  using iterator = super_type::iterator;
  using const_iterator = super_type::const_iterator;
  using value_type = super_type::value_type;

  CouplingConstant() {}
  CouplingConstant(const std::string& expr); 
  CouplingConstant(const strMatrix::row_t& expr_vec);
  ~CouplingConstant() {}
  CouplingConstant& operator=(const std::string& expr); 
  CouplingConstant& operator=(const std::vector<std::string>& expr_vec); 
  CouplingConstant& operator=(const strMatrix& expr_mat); 
  //std::pair<iterator, bool> insert(const value_type& val);
  void create(const unsigned& num_type);
  void add_type(const unsigned& type, const std::string& expr);
  void add_type(const unsigned& type, const std::vector<std::string>& expr_vec);
  void add_type(const unsigned& type, const strMatrix& expr_mat);
  //void create(const value_type& type0, const value_type& type1={0,"_null_"}, 
  //  const value_type& type2={0,"_null_"}, const value_type& type3={0,"_null_"}, 
  //  const value_type& type4={0,"_null_"}, const value_type& type5={0,"_null_"});
  //void add_type(const value_type& val);
  void clear(void); 
  void clear_map(void) { super_type::clear(); } 
  const bool& valid(void) const { return valid_; } 
  const std::string& expression(const unsigned& type) const;
  static const int global_type;
private:
  int num_types_{0}; 
  bool valid_{false};
};

class HamiltonianTerm 
{
public:
  using ccval_t = ComplexMatrix;
  //using BondSiteMap = std::map<unsigned, std::pair<unsigned, unsigned> >;
  HamiltonianTerm() {}
  HamiltonianTerm(const std::string& name, const op::quantum_op& op, const CouplingConstant& cc, 
    const unsigned& size) { construct(name, op, cc, size); }
  ~HamiltonianTerm() {}
  void construct(const std::string& name, const op::quantum_op& op, const CouplingConstant& cc, 
    const unsigned& size);
  void eval_coupling_constant(const ModelParams& cvals, const ModelParams& pvals);
  const op::quantum_op& qn_operator(void) const { return op_; }
  bool is_defined_for(const unsigned& operand_type) const 
    { return is_defined_[operand_type]; }
  const ccval_t& coupling(const unsigned& operand_type) const
    { return cc_values_[operand_type]; }
  const strMatrix& coupling_expr(const unsigned& operand_type) const
    { return cc_.at(operand_type); }
  const std::string& name(void) const { return name_; }
private:
  std::string name_;
  op::quantum_op op_;
  CouplingConstant cc_;
  unsigned max_operand_types_;
  std::vector<bool> is_defined_; // operator defined for a given operand 'type'
  std::vector<ccval_t> cc_values_;
};


} // end namespace model

#endif
