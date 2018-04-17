/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 19:03:43
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-13 22:10:43
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef QUANTUM_OP_H
#define QUANTUM_OP_H

#include <string>

namespace model {

enum class spin { UP, DN, UD, SIGMA, SINGLET };

enum class op_id {
  ni_sigma, ni, cdagc_up, cdagc_dn, cdagc_sigma, sisj, sisj_plus, 
  niup_nidn, cdagup_cdagdn, null
};

enum class op_type { quadratic, pairing, quartic };

namespace op {
class quantum_op 
{
public:
  quantum_op() {}
  quantum_op(const std::string& name, const op_id& id, const spin& s, const op_type& type)
    : name_(name), id_(id), spin_(s), type_(type)
  {
    switch (spin_) {
      case spin::UP:
        spin_up_ = true;
        spin_dn_ = false;
        break;
      case spin::DN:
        spin_up_ = false;
        spin_dn_ = true;
        break;
      case spin::SIGMA:
        spin_up_ = true;
        spin_dn_ = true;
        break;
      default:
        spin_up_ = false;
        spin_dn_ = false;
        break;
    }
    switch (type_) {
      case op_type::quadratic:
        is_quadratic_ = true;
        is_pairing_ = false;
        is_quartic_ = false;
        break;
      case op_type::pairing:
        is_quadratic_ = false;
        is_pairing_ = true;
        is_quartic_ = false;
        break;
      case op_type::quartic:
        is_quadratic_ = false;
        is_pairing_ = false;
        is_quartic_ = true;
        break;
    }
  }
  ~quantum_op() {}
  const std::string& name(void) const { return name_; }
  const op_id& id(void) const { return id_; }
  const spin& sigma(void) const { return spin_; }
  const op_type& type(void) const { return type_; }
  const bool& spin_up(void) const { return spin_up_; }
  const bool& spin_dn(void) const { return spin_dn_; }
  const bool& is_quadratic(void) const { return is_quadratic_; }
  const bool& is_pairing(void) const { return is_pairing_; }
  const bool& is_quartic(void) const { return is_quartic_; }
private:
  std::string name_;
  op_id id_;
  spin spin_;
  op_type type_;
  bool spin_up_;
  bool spin_dn_;
  bool is_quadratic_;
  bool is_pairing_;
  bool is_quartic_;
};

class ni_up : public quantum_op
{
public:
  ni_up() : quantum_op("ni_up", op_id::ni_sigma, spin::UP, op_type::quadratic) {}
};

class ni_dn : public quantum_op
{
public:
  ni_dn() : quantum_op("ni_dn", op_id::ni_sigma, spin::DN, op_type::quadratic) {}
};

// implies both UP or DN
class ni_sigma : public quantum_op
{
public:
  ni_sigma() : quantum_op("ni_sigma", op_id::ni_sigma, spin::SIGMA, op_type::quadratic) {}
};

// implies both UP or DN
class spin_hop : public quantum_op
{
public:
  spin_hop() : quantum_op("spin_hop", op_id::cdagc_sigma, spin::SIGMA, op_type::quadratic) {}
};

class upspin_hop : public quantum_op
{
public:
  upspin_hop() : quantum_op("upspin_hop", op_id::cdagc_up, spin::UP, op_type::quadratic) {}
};

class dnspin_hop : public quantum_op
{
public:
  dnspin_hop() : quantum_op("dnspin_hop", op_id::cdagc_dn, spin::DN, op_type::quadratic) {}
};

class pair_create : public quantum_op
{
public:
  pair_create() : quantum_op("pair_creation", op_id::cdagup_cdagdn, spin::SINGLET, op_type::pairing) {}
};

class sisj_plus : public quantum_op
{
public:
  sisj_plus() : quantum_op("SiSj_ninj/4", op_id::sisj_plus, spin::UD, op_type::quartic) {}
};

class hubbard_int : public quantum_op
{
public:
  hubbard_int() : quantum_op("hubbard", op_id::niup_nidn, spin::UD, op_type::quartic) {}
};

} // end namespace op
} // end namespace model

#endif
