/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-19 11:24:03
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-04-19 11:26:10
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef RBASIS_STATES_H
#define RBASIS_STATES_H

#include <iostream>
#include <string>
#include <cmath>
#include <bitset>
#include <vector>
#include <limits>
#include <stdexcept>

namespace srmf {

class qbitset {
public:
  using size_t = unsigned;
  using long_t = unsigned long long;
  enum {max_bits=std::numeric_limits<long_t>::digits};
  typedef std::bitset<max_bits> bitstring_t;
  qbitset() {}
  qbitset(const unsigned& qval, const unsigned& num_qbits) 
    { construct(qval, num_qbits); }
  ~qbitset() {};
  void construct(const unsigned& qval, const unsigned& num_qbits);
  void reset(void); 
  void reset(const qbitset::size_t& pos); 
  void set(void); 
  void set(const qbitset::size_t& pos); 
  bool raise(const qbitset::size_t& pos); 
  bool lower(const qbitset::size_t& pos); 
  const size_t& size(void) const { return num_bits_; }
  size_t operator[](const unsigned& pos) const;
  size_t bitval(const unsigned& pos) const;
  unsigned long to_ulong() const { return bitstring_.to_ulong(); }
  unsigned long long to_ullong() const { return bitstring_.to_ullong(); }
  friend std::ostream& operator<<(std::ostream& os, const qbitset& qbits_);
private:
  unsigned q_{2};
  unsigned num_bits_{0};
  unsigned block_size_{1};
  std::vector<bitstring_t> bitvals_;
  bitstring_t maxbit_;
  bitstring_t bitstring_;
};

class rotor_basis  
{
public:
  using idx_t = qbitset::long_t;
  using op_result = std::pair<int,idx_t>;
  rotor_basis() {}
  rotor_basis(const unsigned& num_sites, const int& qn_min, const int&qn_max) 
    { construct(num_sites, qn_min, qn_max); }
  ~rotor_basis() {};
  void construct(const unsigned& num_sites, const int& qn_min, const int&qn_max);
  const idx_t& dim(void) const { return ndim_; }
  const idx_t& null_idx(void) const { return null_idx_; }
  op_result apply_cidag(const idx_t& idx, const unsigned& site) const; 
  idx_t apply_cidag(const idx_t& idx) const { return (idx+1)<ndim_ ? idx+1 : null_idx_; }; 
  op_result apply_ni(const idx_t& idx, const unsigned& site) const; 
  int apply_ni(const idx_t& idx) const { return theta_min_ + int(idx); } 
  //int op_ni(const unsigned& site) { return theta_min_; }
  //void random_init(const unsigned int& nsite);
  //void print_out(const int& sfx=-1, std::ostream& os=std::cout) const;
  //unsigned long min_idx(void);
  //unsigned long max_idx(const unsigned& nsite);
  //unsigned long idx(void) const;
private:
  int num_sites_{0};
  int theta_min_{0};
  int theta_max_{0};
  idx_t site_dim_;
  idx_t ndim_;
  idx_t null_idx_;
  std::vector<qbitset::long_t> state_idx_;
  std::vector<qbitset> basis_states_;
};

} // end namespace srmf

#endif