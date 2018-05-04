/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@macbook
* @Date:   2018-04-21 11:41:01
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-05-03 13:56:36
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <cassert>
#include "rbasis_states.h"

namespace srmf {

void qbitset::construct(const unsigned& q, const unsigned& num_qbits)
{
  if (q < 1) {
    throw std::invalid_argument("qbitset::construct: invalid value to argument-1");
  }
  q_ = q;
  num_bits_ = num_qbits;
  if (q_ > 1) {
    unsigned n = q_-1;
    block_size_ = 0;
    while (n) {
      n >>= 1;
      block_size_++;
    }
  }
  //std::cout << "block_size_ = " << block_size_ << "\n";
  if (num_bits_*block_size_ > max_bits) {
    throw std::runtime_error("qbitset::construct: qbitset too long");
  }
  // store the bit values
  bitvals_.clear();
  for (unsigned i=0; i<q_; ++i) { 
    bitvals_.push_back({i}); 
    //std::cout << "bit = " << bitvals_[i] << "\n";
  }
  maxbit_ = bitvals_.back();
  bitstring_.reset();
  //std::cout << bitstring_ << "\n";
}

void qbitset::reset()
{
  bitstring_.reset();
}

void qbitset::reset(const qbitset::size_t& pos)
{
  size_t n = pos * block_size_;
  for (unsigned i=0; i<block_size_; ++i) bitstring_.reset(n+i);
}

void qbitset::set()
{
  for (auto i=0; i<num_bits_; ++i) set(i);
}

void qbitset::set(const qbitset::size_t& pos)
{
  size_t n = pos * block_size_;
  for (unsigned i=0; i<block_size_; ++i) bitstring_[n+i] = maxbit_[i];
}

bool qbitset::raise(const qbitset::size_t& pos) 
{
  assert(pos>=0 || pos<num_bits_);
  size_t n = pos * block_size_;
  bitstring_t bit(0); 
  for (unsigned i=0; i<block_size_; ++i) {
    bit[i] = bitstring_[n+i];
  }
  auto p = bit.to_ulong();
  if (p < bitvals_.size()-1) {
    bit = bitvals_[p+1];
    for (unsigned i=0; i<block_size_; ++i)  bitstring_[n+i] = bit[i];
    //std::cout << "bit = " << bitstring_ << " " << p+1 << "\n";
    return true;
  }
  else return false;
}

bool qbitset::lower(const qbitset::size_t& pos) 
{
  assert(pos>=0 || pos<num_bits_);
  size_t n = pos * block_size_;
  bitstring_t bit(0); 
  for (unsigned i=0; i<block_size_; ++i) {
    bit[i] = bitstring_[n+i];
  }
  auto p = bit.to_ulong();
  if (p > 0) {
    bit = bitvals_[p-1];
    for (unsigned i=0; i<block_size_; ++i) bitstring_[n+i] = bit[i];
    return true;
  }
  else return false;
}

qbitset::size_t qbitset::operator[](const unsigned& pos) const
{
  assert(pos>=0 || pos<num_bits_);
  size_t n = pos * block_size_;
  bitstring_t bit(0); 
  for (unsigned i=0; i<block_size_; ++i) {
    bit[i] = bitstring_[n+i];
  }
  return bit.to_ulong();
}

qbitset::size_t qbitset::bitval(const unsigned& pos) const
{
  assert(pos>=0 || pos<num_bits_);
  size_t n = pos * block_size_;
  bitstring_t bit(0); 
  for (unsigned i=0; i<block_size_; ++i) {
    bit[i] = bitstring_[n+i];
  }
  return bit.to_ulong();
}

std::ostream& operator<<(std::ostream& os, const qbitset& qbits_)
{
  for (int i=qbits_.size()-1; i>=0; --i) {
    size_t n = i * qbits_.block_size_;
    for (int j=qbits_.block_size_-1; j>=0; --j) {
      os << qbits_.bitstring_[n+j];
    }
    os << " ";
  }
  return os;
}


//-------------- rotor_basis --------------
void rotor_basis::construct(const unsigned& num_sites, const int& qn_min, const int&qn_max)
{
  assert(qn_min <= qn_max);
  num_sites_ = num_sites;
  theta_min_ = qn_min;
  theta_max_ = qn_max;
  site_dim_ = (theta_max_-theta_min_)+1;
  ndim_ = std::pow(site_dim_, num_sites_);
  basis_states_.resize(ndim_);
  qbitset qbits_(site_dim_, num_sites_);
  //std::cout << num_sites_ << "\n"; getchar();

  qbits_.set();
  state_idx_.resize(qbits_.to_ullong()+1);
  null_idx_ = ndim_;
  for (auto& idx : state_idx_) idx = null_idx_;

  qbits_.reset();
  basis_states_[0] = qbits_;
  state_idx_[qbits_.to_ullong()] = 0;
  //std::cout << "|" << 0 << "> = " << qbits_ << " " << qbits_.to_ullong() << "\n"; 
  auto max_bit = site_dim_-1;
  for (auto n=1; n<ndim_; ++n) {
    for (auto i=0; i<num_sites_; ++i) {
      if (qbits_.bitval(i) == max_bit) {
        qbits_.reset(i);
      }
      else { 
        qbits_.raise(i);
        break;
      }
    }
    basis_states_[n] = qbits_;
    state_idx_[qbits_.to_ullong()] = n;
    //std::cout << "|" << n << "> = " << qbits_ << " " << qbits_.to_ullong() << "\n"; 
  }
}

rotor_basis::op_result rotor_basis::apply_cidag(const rotor_basis::idx_t& idx, 
  const unsigned& site) const
{
  auto new_state = basis_states_[idx];
  if (new_state.raise(site)) {
    //std::cout << new_state << " " << new_state.to_ullong() << "\n";
    return std::make_pair(1, state_idx_[new_state.to_ullong()]);
  }
  else return std::make_pair(1, null_idx_);
}

rotor_basis::op_result rotor_basis::apply_ni(const rotor_basis::idx_t& idx, 
  const unsigned& site) const
{
  int ntheta = theta_min_ + basis_states_[idx].bitval(site);
  return std::make_pair(ntheta, idx);
}

} // end namespace srmf
