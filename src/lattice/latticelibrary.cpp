/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2016-01-17 21:32:15
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-17 11:50:38
*----------------------------------------------------------------------------*/
#include <stdexcept>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include "lattice.h"
//#include "graph.h"

namespace lattice {

// define the lattice
int Lattice::define_lattice(void) 
{
  using pos = Eigen::Vector3i;
  using vec = Eigen::Vector3d;
  unsigned type, src, tgt;
  vec a1, a2, a3, coord;
  pos src_offset, tgt_offset, cell;

  /*------------- 'SQUARE' lattice--------------*/
  if (lname == "SQUARE") {
    /*
    // type
    lid = lattice_id::SQUARE;
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open};
    // basis vectors
    set_basis_vectors(a1=vec(1,0,0), a2=vec(0,1,0), a3=vec(0,0,0));
    // add sites
    add_basis_site(type=0, coord=vec(0,0,0));
    // add bonds
    add_bond(type=0, ngb=1, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(1,0,0));
    add_bond(type=1, ngb=1, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(0,1,0));
    */
  }

  else if (lname == "SQUARE_2BAND") {
    // type
    lid = lattice_id::SQUARE_2BAND;
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open};
    // basis vectors
    set_basis_vectors(a1=vec(1,0,0), a2=vec(0,1,0), a3=vec(0,0,0));
    // add sites
    int orbitals;
    add_basis_site(orbitals=2, coord=vec(0,0,0));
    add_bond(type=0, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(1,0,0));
    add_bond(type=1, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(0,1,0));

    // bonds along +x
    /*
    // add sites
    add_basis_site(type=0, atom=0, coord=vec(0,0,0));
    add_basis_site(type=1, atom=0, coord=vec(0,0,0));
    add_bond(type=0, ngb=1, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(1,0,0));
    add_bond(type=1, ngb=1, src=0, src_offset=pos(0,0,0), tgt=1, tgt_offset=pos(1,0,0));
    add_bond(type=2, ngb=1, src=1, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(1,0,0));
    add_bond(type=3, ngb=1, src=1, src_offset=pos(0,0,0), tgt=1, tgt_offset=pos(1,0,0));
    // bonds along +y
    add_bond(type=0, ngb=1, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(0,1,0));
    add_bond(type=4, ngb=1, src=0, src_offset=pos(0,0,0), tgt=1, tgt_offset=pos(0,1,0));
    add_bond(type=5, ngb=1, src=1, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(0,1,0));
    add_bond(type=3, ngb=1, src=1, src_offset=pos(0,0,0), tgt=1, tgt_offset=pos(0,1,0));
    */
  }

  else if (lname == "SIMPLE_CUBIC") {
    // type
    /*
    lid = lattice_id::SIMPLECUBIC;
    // basis vectors
    set_basis_vectors(a1=vec(1,0,0), a2=vec(0,1,0), a3=vec(0,0,1));
    // add sites
    add_basis_site(type=0, coord=vec(0,0,0));
    // add bonds
    add_bond(type=0, ngb=1, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(1,0,0));
    add_bond(type=0, ngb=1, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(0,1,0));
    add_bond(type=0, ngb=1, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(0,0,1));
    */
  }

  /*------------- 'CHAIN' lattice--------------*/
  else if (lname == "CHAIN") {
    /*
    lid = lattice_id::CHAIN;
    extent[dim2] = Extent{1, boundary_type::open, boundary_type::open};
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open};
    // add sites
    add_basis_site(type=0, coord=vec(0,0,0));
    // add bonds
    add_bond(type=0, ngb=1, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(1,0,0));
    */
  }

  else if (lname == "HONEYCOMB") {
    /*
    // type
    lid = lattice_id::HONEYCOMB;
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open};
    // basis vectors
    double x = 0.5; 
    double y = -0.5*std::sqrt(3.0); 
    set_basis_vectors(a1=vec(x,y,0), a2=vec(x,-y,0), a3=vec(0,0,0));
    // sites
    add_basis_site(type=0, coord=vec(0,0,0));
    y = 1.0/std::sqrt(3.0); 
    add_basis_site(type=1, coord=vec(0,y,0));
    // NN bonds
    add_bond(type=0, ngb=1, src=0, src_offset=pos(0,0,0), tgt=1, tgt_offset=pos(0,0,0));
    add_bond(type=0, ngb=1, src=0, src_offset=pos(0,0,0), tgt=1, tgt_offset=pos(1,0,0));
    add_bond(type=0, ngb=1, src=1, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(0,1,0));
    // NNN bonds
    add_bond(type=1, ngb=2, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(1,0,0));
    add_bond(type=1, ngb=2, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(0,1,0));
    add_bond(type=2, ngb=2, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(1,1,0));
    add_bond(type=3, ngb=2, src=1, src_offset=pos(0,0,0), tgt=1, tgt_offset=pos(1,0,0));
    add_bond(type=3, ngb=2, src=1, src_offset=pos(0,0,0), tgt=1, tgt_offset=pos(0,1,0));
    add_bond(type=4, ngb=2, src=1, src_offset=pos(0,0,0), tgt=1, tgt_offset=pos(1,1,0));
    */
  }

  else if (lname == "KAGOME") {
    /*
    // type
    lid = lattice_id::KAGOME;
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open};
    // basis vectors
    double x = 1.0; 
    double y = std::sqrt(3.0); 
    set_basis_vectors(a1=vec(x,y,0), a2=vec(-x,y,0), a3=vec(0,0,0));
    // add sites
    // A-site
    add_basis_site(type=0, coord=vec(0,0,0));
    // B-site
    x = 0.5; y = 0.5*std::sqrt(3.0);
    add_basis_site(type=1, coord=vec(x,y,0));
    // C-site
    add_basis_site(type=2, coord=vec(-x,y,0));
    // NN-1 bonds
    // A-B (intra-cell) 
    add_bond(type=0,ngb=1,src=0,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,0,0));
    // A-C (intra-cell) 
    add_bond(type=1,ngb=1,src=0,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(0,0,0));
    // B-C (intra-cell) 
    add_bond(type=0,ngb=1,src=1,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(0,0,0));
    // B-A (inter-cell) 
    add_bond(type=1,ngb=1,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    // C-A (inter-cell) 
    add_bond(type=0,ngb=1,src=2,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));
    // B-C (inter-cell) 
    add_bond(type=0,ngb=1,src=1,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(1,-1,0));

    // NN-2 bonds
    // B-C (intra-cell) 
    add_bond(type=3,ngb=2,src=1,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(1,0,0));
    // C-A (inter-cell) 
    add_bond(type=3,ngb=2,src=2,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    // B-A (inter-cell) 
    add_bond(type=2,ngb=2,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));
    // C-B (inter-cell) 
    add_bond(type=2,ngb=2,src=2,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,1,0));
    // A-C (intra-cell) 
    add_bond(type=2,ngb=2,src=0,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(1,-1,0));
    // B-A (inter-cell) 
    add_bond(type=2,ngb=2,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,-1,0));
    */
  }

  else if (lname == "PYROCHLORE_V1") {
    /*
    // type
    lid = lattice_id::PYROCHLORE_V1;
    // basis vectors
    double x = 1.0; 
    double y = std::sqrt(3.0); 
    double y2 = 2.0/std::sqrt(3.0);
    double z = std::sqrt(8.0/3.0);
    set_basis_vectors(a1=vec(x,y,0), a2=vec(-x,y,0), a3=vec(0,y2,z));
    // add sites
    // A-site
    add_basis_site(type=0, coord=vec(0,0,0));
    // B-site
    x = 0.5; y = 0.5*std::sqrt(3.0);
    add_basis_site(type=1, coord=vec(x,y,0));
    // C-site
    add_basis_site(type=2, coord=vec(-x,y,0));
    // D-site
    y = 1.0/std::sqrt(3.0);
    z = std::sqrt(2.0/3.0);
    add_basis_site(type=3, coord=vec(0,y,z));

    // NN-1 bonds
    // A-B (intra-cell) 
    add_bond(type=0,ngb=1,src=0,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,0,0));
    // A-C (intra-cell) 
    add_bond(type=1,ngb=1,src=0,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(0,0,0));
    // B-C (intra-cell) 
    add_bond(type=0,ngb=1,src=1,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(0,0,0));
    // B-A (inter-cell) 
    add_bond(type=1,ngb=1,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    // C-A (inter-cell) 
    add_bond(type=0,ngb=1,src=2,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));
    // B-C (inter-cell) 
    add_bond(type=0,ngb=1,src=1,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(1,-1,0));

    // NN-2 bonds
    // B-C (intra-cell) 
    add_bond(type=3,ngb=2,src=1,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(1,0,0));
    // C-A (inter-cell) 
    add_bond(type=3,ngb=2,src=2,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    // B-A (inter-cell) 
    add_bond(type=2,ngb=2,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));
    // C-B (inter-cell) 
    add_bond(type=2,ngb=2,src=2,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,1,0));
    // A-C (intra-cell) 
    add_bond(type=2,ngb=2,src=0,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(1,-1,0));
    // B-A (inter-cell) 
    add_bond(type=2,ngb=2,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,-1,0));

    // z-dir bonds
    // A-D (intra-cell) 
    add_bond(type=4,ngb=1,src=0,src_offset=pos(0,0,0),tgt=3,tgt_offset=pos(0,0,0));
    // B-D (intra-cell) 
    add_bond(type=4,ngb=1,src=1,src_offset=pos(0,0,0),tgt=3,tgt_offset=pos(0,0,0));
    // C-D (intra-cell) 
    add_bond(type=4,ngb=1,src=2,src_offset=pos(0,0,0),tgt=3,tgt_offset=pos(0,0,0));
    // D-A (inter-cell) 
    add_bond(type=4,ngb=1,src=3,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,0,1));
    // D-B (inter-cell) 
    add_bond(type=4,ngb=1,src=3,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(-1,0,1));
    // D-C (inter-cell) 
    add_bond(type=4,ngb=1,src=3,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(0,-1,1));
    */
  }

  else if (lname == "PYROCHLORE_3D") {
    // type
    lid = lattice_id::PYROCHLORE_3D;
    // basis vectors
    set_basis_vectors(a1=vec(0,2,2), a2=vec(2,0,2), a3=vec(2,2,0));
    //set_basis_vectors(a1=vec(2,0,2), a2=vec(0,2,2), a3=vec(2,2,0));
    // add sites
    int orbitals;
    // Ir site-0
    add_basis_site(type=0, orbitals=6, coord=vec(0,0,0));
    add_basis_site(type=1, orbitals=6, coord=0.5*basis_vector_a1());
    add_basis_site(type=2, orbitals=6, coord=0.5*basis_vector_a2());
    add_basis_site(type=3, orbitals=6, coord=0.5*basis_vector_a3());

    // Intra-cell bonds
    add_bond(type=0,src=0,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,0,0));
    add_bond(type=1,src=0,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(0,0,0));
    add_bond(type=2,src=0,src_offset=pos(0,0,0),tgt=3,tgt_offset=pos(0,0,0));
    add_bond(type=3,src=1,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(0,0,0));
    add_bond(type=4,src=1,src_offset=pos(0,0,0),tgt=3,tgt_offset=pos(0,0,0));
    add_bond(type=5,src=2,src_offset=pos(0,0,0),tgt=3,tgt_offset=pos(0,0,0));

    // Inter-cell bonds
    add_bond(type=6,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    add_bond(type=7,src=2,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));
    add_bond(type=8,src=3,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,0,1));

    add_bond(type=3,src=1,src_offset=pos(0,0,0),tgt=2,tgt_offset=pos(1,-1,0));
    add_bond(type=4,src=1,src_offset=pos(0,0,0),tgt=3,tgt_offset=pos(1,0,-1));
    add_bond(type=5,src=2,src_offset=pos(0,0,0),tgt=3,tgt_offset=pos(0,1,-1));
  }


  /*------------- undefined lattice--------------*/
  else {
    throw std::range_error("error: latticelibrary: undefined lattice");
  }
  return 0;
}

// read lattice parameters
int Lattice::construct(const input::Parameters& parms) 
{

  int info;
  // name
  lname = parms.set_value("lattice", "NULL");
  boost::to_upper(lname);

  // sizes
  for (unsigned dim=dim1; dim<=dim3; ++dim) {
    std::string lsize = "lsize" + std::to_string(dim+1);
    extent[dim].size = parms.set_value(lsize, 1, info);
    if (extent[dim].size<1) throw std::range_error("error: latticelibrary: invalid lattice size");
  }

  // boundary conditions
  std::string bc; 
  for (unsigned dim=dim1; dim<=dim3; ++dim) {
    std::string lbc = "bc" + std::to_string(dim+1);
    bc = parms.set_value(lbc, "open", info);
    extent[dim].periodicity = boundary_condition(bc);
    extent[dim].bc = extent[dim].periodicity;
    if (extent[dim].bc == boundary_type::antiperiodic) extent[dim].bc = boundary_type::periodic;
  }

  // empty unitcell
  Unitcell::clear();

  // impurities
  //impurity_sites_.clear();
  //impurity_bonds_.clear();

  define_lattice();
  finalize_lattice();
  //construct_graph();

  return 0;
}

int Lattice::finalize_lattice(void) 
{
  // Finalize the unit cell
  Unitcell::finalize();

  // copy the user set dimensions
  for (unsigned dim=dim1; dim<=dim3; ++dim) copy_extent[dim] = extent[dim];

  // is it necessary to construct 'symmetrized lattice'
  for (unsigned dim=dim1; dim<=dim3; ++dim) {
    if (extent[dim].bc==boundary_type::open && extent[dim].size>1) {
      symmetrize_lattice();
      break;
    }
  }

  // number of unit cells & sites
  num_layer_cells_ = extent[dim1].size * extent[dim2].size;
  num_total_cells_ = num_layer_cells_ * extent[dim3].size;
  num_basis_sites_ = Unitcell::num_sites();
  num_total_sites_ = num_total_cells_ * num_basis_sites_;

  // indexing of the 'basis orbitals'
  basis_index_.clear();
  basis_index_.resize(num_basis_sites_);
  unsigned idx = 0;
  for (unsigned i=0; i<num_basis_sites_; ++i) {
    for (unsigned j=0; j<basis_site(i).num_orbitals(); ++j) {
      basis_index_[i].push_back(idx++);
    }
  }

  // check
  /*std::cout << "------Sites-------\n";
  for (unsigned i=0; i<unitcell.num_site(); ++i) {
    std::cout << Unitcell::site(i) << std::endl;
  }
  std::cout << "------Bonds-------\n";
  for (unsigned i=0; i<unitcell.num_bond(); ++i) {
    std::cout << Unitcell::bond(i) << std::endl;
  }*/

  // 'vector' & 'vector_id' attributes of the bonds
  std::map<int,unsigned> vecid_map;
  unsigned id=0;
  for (unsigned i=0; i<Unitcell::num_bonds(); ++i) {
    Vector3i ivec = Unitcell::bond(i).tgt().bravindex()-Unitcell::bond(i).src().bravindex();
    int key = ivec[0] + ivec[1]*extent[dim1].size + ivec[2]*num_layer_cells_;
    auto it = vecid_map.find(key);
    if (it != vecid_map.end()) Unitcell::bond(i).set_vector_id(it->second);
    else {
      vecid_map.insert({key, id});
      Unitcell::bond(i).set_vector_id(id);
      id++;
    }
    Unitcell::bond(i).set_vector((ivec[0]*vector_a1()+ivec[1]*vector_a2()+ivec[2]*vector_a3()));
    //std::cout << "bond " << i << ": vector_id = " << bond(i).vector_id() << "\n";
  }

  return 0;
}

int Lattice::symmetrize_lattice(void) 
{
  // initially, the 'dim' with periodic bc has size = 1
  spatial_dim = 0;
  Vector3d bvec;
  for (unsigned dim=dim1; dim<=dim3; ++dim) {
    if (extent[dim].bc==boundary_type::periodic) {
      spatial_dim++;
      // temporarily set size = 1 for dim with PBC
      extent[dim].size = 1;
      switch (dim) {
        case dim1: bvec = basis_vector_a1(); break;
        case dim2: bvec = basis_vector_a2(); break;
        case dim3: bvec = basis_vector_a3(); break;
      }
    }
  }
  // if 1 dimensional lattice, rotate the lattice to align 'bvec' along x-direction
  if (spatial_dim == 1) {
    // rotation matrix to do that
    Eigen::Matrix3d matrix = rotation_matrix(bvec, Vector3d(1.0, 0.0, 0.0));
    // rotate the unitcell
    Unitcell::rotate_by(matrix);
  }

  // number of unit cells & sites
  num_layer_cells_ = extent[dim1].size * extent[dim2].size;
  num_total_cells_ = num_layer_cells_ * extent[dim3].size;
  num_basis_sites_ = Unitcell::num_sites();
  num_total_sites_ = num_total_cells_ * num_basis_sites_;

  // Add the sites & the bonds to the symmetrized unitcell
  std::vector<Site> sites;
  std::vector<Bond> bonds;
  Unitcell translated_cell;
  Vector3i bravindex(0,0,0);
  for (unsigned i=0; i<num_total_cells_; ++i) {
    translated_cell = get_translated_cell(bravindex);
    // collect the sites
    for (unsigned n=0; n<translated_cell.num_sites(); ++n) 
      sites.push_back(translated_cell.site(n));
    // collect the bonds
    for (unsigned n=0; n<translated_cell.num_bonds(); ++n) 
      bonds.push_back(translated_cell.bond(n));
    bravindex = get_next_bravindex(bravindex);
  }

  // replace the old sites & bonds
  Unitcell::clear_sites();
  unsigned i = 0;
  for (auto& s : sites) {
    s.reset_uid(i++); 
    s.reset_bravindex(Vector3i(0,0,0));
    s.reset_cell_coord(Vector3d(0.0,0.0,0.0));
    Unitcell::add_site(s);
    //std::cout<<"site = "<< s <<"\n"; getchar();
  }
  Unitcell::clear_bonds();
  for (auto& b : bonds) {
    //std::cout<<"bond = "<<b.src()<<" "<<b.tgt()<<"\n";
    //std::cout<<"bond = "<<b.bravindex().transpose()<<"\n";
    if (connect_bond(b, sites)) {
      b.reset_bravindex(Vector3i(0,0,0));
      Unitcell::add_bond(b);
      //std::cout << "***connected*****\n";
      //std::cout<<"bond = "<<b.src()<<" "<<b.tgt()<<"\n"; getchar();
    }
    ///else std::cout << "***NOT connected****\n";
    ///getchar();
  }


  //std::cout << unitcell.vector_a1() << "\n";
  //std::cout << unitcell.vector_a2() << "\n";
  //std::cout << unitcell.vector_a3() << "\n";

  // extent & basis vectors of the symmetrized lattice
  for (unsigned dim=dim1; dim<=dim3; ++dim) {
    extent[dim] = copy_extent[dim];
    if (extent[dim].bc != boundary_type::periodic) {
      extent[dim].size = 1;
      switch (dim) {
        case dim1: Unitcell::reset_a1(Vector3d(0.0,0.0,0.0)); break;
        case dim2: Unitcell::reset_a2(Vector3d(0.0,0.0,0.0)); break;
        case dim3: Unitcell::reset_a3(Vector3d(0.0,0.0,0.0)); break;
      }
    }
  }
  //std::cout << unitcell.vector_a1() << "\n";
  //std::cout << unitcell.vector_a2() << "\n";
  //std::cout << unitcell.vector_a3() << "\n";

  return 0;
}






} // end namespace lattice
