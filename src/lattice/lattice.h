/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-01-25 18:05:03
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-17 11:41:52
*----------------------------------------------------------------------------*/
#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <string>
#include <tuple>
#include <set>
#include <stdexcept>
#include "../scheduler/task.h"
#include <Eigen/Dense>
#include <boost/graph/adjacency_list.hpp>
#include "constants.h"

namespace lattice {

// Global Constants
const unsigned MAX_SITE_TYPES = 20;
const unsigned MAX_BOND_TYPES = 40;

/*---------------lattice types-----------------*/
enum class lattice_id {
  UNDEFINED, SQUARE, SQUARE_2BAND, CHAIN, HONEYCOMB, SIMPLECUBIC, 
  KAGOME, PYROCHLORE_V1, PYROCHLORE_3D
};

/*---------------Lattice site class-----------------*/
using Vector3i = Eigen::Vector3i;
using Vector3d = Eigen::Vector3d;
//using BrvaisIdx = Eigen::Matrix<unsigned, 3, 1>;

class Site 
{
public:
  // ctors
  Site() {}
  Site(const unsigned& uid, const unsigned& type, const unsigned& atomid, const Vector3i& bravindex, 
    const Vector3d& coord, const Vector3d& cell_coord);
  ~Site() {}
  // setter functions
  static void reset_count(void) { num_site=0; }
  void reset_type(const unsigned& t) { type_=t; }
  void reset_uid(const unsigned& uid) { uid_=uid; }
  void reset_bravindex(const Vector3i& idx) { bravindex_=idx; }
  void reset_coord(const Vector3d& v) { coord_=v; }
  void reset_cell_coord(const Vector3d& v) { cell_coord_=v; }
  void translate_by(const int& id_offset, const Vector3i& bravindex_offset, const Vector3d& coord_offset); 

  // getter functions
  const int& id(void) const {return id_;}
  const unsigned& uid(void) const {return uid_;}
  const unsigned& type(void) const {return type_;}
  const unsigned& atomid(void) const {return atomid_;}
  const unsigned& num_orbitals(void) const {return num_orbitals_;}
  const Vector3i& bravindex(void) const { return bravindex_; }
  const Vector3d& coord(void) const {return coord_;}
  const Vector3d& cell_coord(void) const {return cell_coord_;}
  // friends
  friend std::ostream& operator<<(std::ostream& os, const Site& site);
private:
  static unsigned num_site;
  int id_ {0};
  unsigned uid_ {0}; // local id within a unitcell
  unsigned type_ {0};
  unsigned atomid_ {0};
  unsigned num_orbitals_ {1};
  Vector3i bravindex_ {Vector3i(0, 0, 0)};
  Vector3d coord_ {Vector3d(0.0, 0.0, 0.0)};
  Vector3d cell_coord_ {Vector3d(0.0, 0.0, 0.0)};
};

/*---------------Lattice bond class-----------------*/
class Bond : public std::pair<Site, Site> 
{
public:
  // ctors
  Bond() {}
  Bond(const unsigned& type, const unsigned& ngb, const Site& src, const Vector3i& bravindex,
    const Site& tgt, const Vector3i& tgt_offset);
  //Bond(const unsigned& type, const unsigned& ngb, const Vector3i& bravindex, const unsigned& src_id, 
  //  const Vector3i& src_offset, const unsigned& tgt_id, const Vector3i& tgt_offset, const int& sign);
  ~Bond() {}
  // setter functions
  static void reset_count(void) { num_bond=0; }
  void reset_id(const unsigned& id) { id_=id; }
  void reset_type(const unsigned& t) { type_=t; }
  //void reset_src_offset(const Vector3i& idx) { src_offset_=idx; }
  //void reset_tgt_offset(const Vector3i& idx) { tgt_offset_=idx; }
  void reset_bravindex(const Vector3i& idx) { bravindex_=idx; }
  void set_vector_id(const unsigned& id) { vector_id_=id; }
  void set_vector(const Vector3d& R) { vector_=R; }
  //void shift_target_ids(const int& id_offset) { src_ += id_offset; tgt_ += id_offset; }
  void translate_by(const Vector3i& bravindex_offset) { bravindex_ += bravindex_offset; } 
  void connect(const Site& src, const Vector3i& src_offset, const Site& tgt, 
    const Vector3i& tgt_offset, const int& sign);
  //void connect(const unsigned& src_id, const Vector3i& src_offset, const unsigned& tgt_id, 
  //  const Vector3i& tgt_offset, const int& sign);
  // getter functions
  const int& id(void) const { return id_; }
  const unsigned& type(void) const {return type_;}
  unsigned src_id(void) const { return first.id(); }
  unsigned tgt_id(void) const { return second.id(); }
  const Site& src(void) const { return first; }
  const Site& tgt(void) const { return second; }
  const unsigned& vector_id(void) const { return vector_id_; }
  const Vector3d& vector(void) const { return vector_; }
  int sign(void) const { return sign_; }
  Vector3i bravindex(void) const { return bravindex_; }
  //Vector3i src_offset(void) const { return src_offset_; }
  //Vector3i tgt_offset(void) const { return tgt_offset_; }
  // friends
  friend std::ostream& operator<<(std::ostream& os, const Bond& bond);
private:
  static unsigned num_bond;
  int id_ {0};
  unsigned type_ {0};
  unsigned ngb_ {0};
  unsigned vector_id_{0}; // integer id for the following vector
  Vector3d vector_{Vector3d(0,0,0)}; // coordinate of 'tgt cell' wrt 'src cell'
  //unsigned src_ {0}; 
  //unsigned tgt_ {0}; 
  int sign_ {1}; // = -1 if across an antiperiodic boundary
  Vector3i bravindex_ {Vector3i(0, 0, 0)};
  //Vector3i src_offset_ {Vector3i(0, 0, 0)};
  //Vector3i tgt_offset_ {Vector3i(0, 0, 0)};
};

/*---------------Unitcell class-----------------*/
class Unitcell 
{
public:
  // ctors
  Unitcell() {}
  ~Unitcell() {}
  // setter functions
  int add_site(const unsigned& type, const unsigned& num_orbitals, const Vector3d& site_coord); 
  int add_site(const Site& s) { sites.push_back(s); return sites.back().id(); }
  int add_site(const Site& s, const Vector3i& bravindex, const Vector3d& cell_coord);
  int add_bond(const Bond& b) { bonds.push_back(b); return bonds.back().id(); }
  int add_bond(const unsigned& type, const unsigned& ngb, const unsigned& src_id, const Vector3i& src_offset,
    const unsigned& tgt_id, const Vector3i& tgt_offset); 
  void set_basis(const Vector3d& av1, const Vector3d& av2, const Vector3d& av3);
  void reset_a1(const Vector3d& av) { a1=av; }
  void reset_a2(const Vector3d& av) { a2=av; }
  void reset_a3(const Vector3d& av) { a3=av; }
  void clear(void); 
  void clear_sites(void) { sites.clear(); }
  void clear_bonds(void) { bonds.clear(); }
  void finalize(void);
  // getter functions
  unsigned num_sites(void) const { return sites.size(); }
  unsigned num_bonds(void) const { return bonds.size(); }
  const unsigned& num_orbitals(void) const { return num_orbitals_; }
  const Vector3d& vector_a1(void) const { return a1; }
  const Vector3d& vector_a2(void) const { return a2; }
  const Vector3d& vector_a3(void) const { return a3; }
  const Vector3i& bravindex(void) const { return bravindex_; }
  const Vector3d& coord(void) const {return coord_;}
  const Site& site(const unsigned& i) const { return sites[i]; }
  const Bond& bond(const unsigned& i) const { return bonds[i]; }
  Bond& bond(const unsigned& i) { return bonds[i]; }
  unsigned num_site_types(void) const { return sitetypes_map_.size(); }
  unsigned num_bond_types(void) const { return bondtypes_map_.size(); }
  std::map<unsigned,unsigned> sitetypes_map(void) const { return sitetypes_map_; }
  std::map<unsigned,unsigned> bondtypes_map(void) const { return bondtypes_map_; }
  // transform
  void translate_by(const Vector3i& bravindex_offset, const int& cell_id_offset);
  void rotate_by(const Eigen::Matrix3d& matrix);

private:
  int id {0};
  unsigned max_site_type_val {0};
  unsigned max_bond_type_val {0};
  unsigned max_neighb_val {0};
  Vector3d a1 {Vector3d(0.0, 0.0, 0.0)};
  Vector3d a2 {Vector3d(0.0, 0.0, 0.0)};
  Vector3d a3 {Vector3d(0.0, 0.0, 0.0)};
  Vector3i bravindex_ {Vector3i(0, 0, 0)};
  Vector3d coord_ {Vector3d(0.0, 0.0, 0.0)};
  std::vector<Site> sites;
  std::vector<Bond> bonds;
  unsigned num_orbitals_ {0};
  std::map<unsigned,unsigned> sitetypes_map_; // user set value to contiguous value 
  std::map<unsigned,unsigned> bondtypes_map_; // user set value to contiguous value 
};

/*---------------spatial dimension type-----------------*/
enum class boundary_type {open, periodic, antiperiodic};

class Lattice : public Unitcell
{
public:
  // ctors
  Lattice() {}
  Lattice(const input::Parameters& parms) { construct(parms); }
  ~Lattice() {}
  // setter functions
  int construct(const input::Parameters& parms);
  void set_basis_vectors(const Vector3d& av1, const Vector3d& av2, const Vector3d& av3);
  int add_basis_site(const Vector3d& site_coord);
  int add_basis_site(const unsigned& num_orbitals, const Vector3d& site_coord);
  int add_basis_site(const unsigned& type, const unsigned& num_orbitals, 
    const Vector3d& site_coord);
  int add_bond(const unsigned& src_id, const Vector3i& src_offset,
    const unsigned& tgt_id, const Vector3i& tgt_offset); 
  int add_bond(const unsigned& type, const unsigned& src_id, const Vector3i& src_offset,
    const unsigned& tgt_id, const Vector3i& tgt_offset); 
  int add_bond(const unsigned& type, const unsigned& ngb, const unsigned& src_id, const Vector3i& src_offset,
    const unsigned& tgt_id, const Vector3i& tgt_offset); 
  int add_bond(const Bond& b); 

  // getter functions
  std::string name(void) const { return lname; }
  lattice_id id(void) const { return lid; }
  unsigned num_basis_sites(void) const { return Unitcell::num_sites(); }
  unsigned num_basis_bonds(void) const { return Unitcell::num_bonds(); }
  unsigned num_basis_orbitals(void) const { return Unitcell::num_orbitals(); }
  const unsigned& dimension(void) const { return spatial_dim; }
  const unsigned& num_sites(void) const { return num_total_sites_; }
  const unsigned& num_unitcells(void) const { return num_total_cells_; }
  Vector3d basis_vector_a1(void) const { return Unitcell::vector_a1(); }
  Vector3d basis_vector_a2(void) const { return Unitcell::vector_a2(); }
  Vector3d basis_vector_a3(void) const { return Unitcell::vector_a3(); }
  int size1(void) const { return static_cast<int>(extent[dim1].size); }
  int size2(void) const { return static_cast<int>(extent[dim2].size); }
  int size3(void) const { return static_cast<int>(extent[dim3].size); }
  boundary_type bc1(void) const { return extent[dim1].bc; }
  boundary_type bc2(void) const { return extent[dim2].bc; }
  boundary_type bc3(void) const { return extent[dim3].bc; }
  //boundary_type bc2(void) const { return extent[dim2].bc; }
  //boundary_type bc3(void) const { return extent[dim3].bc; }
  boundary_type bc1_periodicity(void) const { return extent[dim1].periodicity; }
  boundary_type bc2_periodicity(void) const { return extent[dim2].periodicity; }
  boundary_type bc3_periodicity(void) const { return extent[dim3].periodicity; }

  // other methods 
  //Vector3i boundary_wrap(const Vector3i& cell_idx) const;
  std::pair<Vector3i, int> boundary_wrap(const Vector3i& cell_idx) const;
  Vector3i get_next_bravindex(const Vector3i& current_index) const;
  Unitcell get_translated_cell(const Vector3i& bravindex_offset) const;
  int mapped_site_id(const unsigned& local_id, const Vector3i& bravindex) const;
  unsigned translation_mapped_site(const unsigned& uid, const Vector3i& bravindex,
    const Vector3i& translation_vec) const;
  bool connect_bond(Bond& bond, const std::vector<Site>& sites) const;
  const Site& basis_site(const unsigned& i) const { return Unitcell::site(i); }
  const Bond& basis_bond(const unsigned& i) const { return Unitcell::bond(i); }
  const unsigned& basis_index_number(const unsigned& site, const unsigned& orbital) const
    { return basis_index_[site][orbital]; }

  // for lattices with doped impurities
  //void add_new_bondtype(const unsigned& type);
private:
  struct Extent {unsigned size; boundary_type bc; boundary_type periodicity;};
  enum Dimension {dim1, dim2, dim3};
  lattice_id lid {lattice_id::SQUARE};
  std::string lname {""};
  unsigned spatial_dim {0};

  // lattice dimensions
  Extent extent[3] = {Extent{1, boundary_type::open, boundary_type::open}, 
                      Extent{1, boundary_type::open, boundary_type::open},
                      Extent{1, boundary_type::open, boundary_type::open}
                     };

  // copy of user-set lattice dimensions
  Extent copy_extent[3] {Extent{1, boundary_type::open, boundary_type::open}, 
                         Extent{1, boundary_type::open, boundary_type::open},
                         Extent{1, boundary_type::open, boundary_type::open}
                        };
  
  // number of unit cells in total and in one layer (for symmetrized lattice)
  unsigned num_total_cells_{1};
  unsigned num_layer_cells_{1};
  unsigned num_total_sites_{0};
  unsigned num_basis_sites_{0};

  // basis state index numbering
  std::vector<std::vector<unsigned> > basis_index_;

  // for lattices with impurities
  std::vector<Site> impurity_sites_;
  std::vector<Bond> impurity_bonds_;

  // helper functions
  int define_lattice(void); 
  int finalize_lattice(void); 
  int symmetrize_lattice(void);
  //int construct_graph(void); 
  boundary_type boundary_condition(std::string& bc) const;
  Eigen::Matrix3d rotation_matrix(const Eigen::Vector3d& r, const Eigen::Vector3d& r_prime);
};


} // end namespace lattice

#endif







