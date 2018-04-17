/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-01 21:13:27
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-16 00:17:18
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "blochbasis.h"

namespace basis {

BlochBasis::BlochBasis(const lattice::LatticeGraph& graph) 
{
  construct(graph);
}

int BlochBasis::construct(const lattice::LatticeGraph& graph)
{
  /*
  if (disorded_system) {
    // no translational symmetry
    // only k=0 point
    num_kpoint_ = 1;
    this->clear();
    this->push_back(Vector3d(0.0,0.0,0.0));
    // subspace basis
    subspace_dimension_ = graph.num_sites();
    subspace_basis_.resize(subspace_dimension_);
    for (unsigned i=0; i<subspace_dimension_; ++i) {
      basis_state s = graph.site(i);
      subspace_basis_[i] = s; 
    }
    null_idx_ = subspace_basis_.size();
    // index of the 'representative state' of a site
    representative_state_idx_.resize(graph.num_sites());
    for (unsigned i=0; i<graph.num_sites(); ++i) {
      representative_state_idx_[i] = i;
    }
    return 0;
  }
  */
  make_kpoints(graph.lattice());
  make_subspace_basis(graph);
  return 0;
}

void BlochBasis::make_kpoints(const lattice::Lattice& lattice)
{
  Vector3d a1, a2, a3;
  double v;
  using bc = lattice::boundary_type;

  // reciprocal lattice vectors 
  b1 = Vector3d(0.0,0.0,0.0);
  b2 = Vector3d(0.0,0.0,0.0);
  b3 = Vector3d(0.0,0.0,0.0);

  unsigned symmetry_type = 0;
  if (lattice.bc1() == bc::periodic) {
    a1 = lattice.vector_a1();
    b1 = two_pi() * a1 / a1.dot(a1); 
    symmetry_type = symmetry_type + 1;
  } 

  if (lattice.bc2() == bc::periodic) {
    a1 = lattice.vector_a2();
    a2 = lattice.vector_a2();
    switch (symmetry_type) {
      case 0:
        b2 = two_pi() * a2 / a2.dot(a2); 
        break;
      case 1:
        a3 = a1.cross(a2);
        v = a1.dot(a2.cross(a3));
        b1 = two_pi() * a2.cross(a3) / v;
        b2 = two_pi() * a3.cross(a1) / v;
        break;
      default: break;
    }
    symmetry_type = symmetry_type + 2;
  } 

  if (lattice.bc3() == bc::periodic) {
    a1 = lattice.vector_a1();
    a2 = lattice.vector_a2();
    a3 = lattice.vector_a3();
    switch (symmetry_type) {
      case 0:
        b3 = two_pi() * a3 / a3.dot(a3); 
        break;
      case 1:
        a2 = a3.cross(a1);
        v = a1.dot(a2.cross(a3));
        b1 = two_pi() * a2.cross(a3) / v;
        b3 = two_pi() * a1.cross(a2) / v;
        break;
      case 2:
        a1 = a2.cross(a3);
        v = a1.dot(a2.cross(a3));
        b2 = two_pi() * a3.cross(a1) / v;
        b3 = two_pi() * a1.cross(a2) / v;
        break;
      case 3:
        v = a1.dot(a2.cross(a3));
        b1 = two_pi() * a2.cross(a3) / v;
        b2 = two_pi() * a3.cross(a1) / v;
        b3 = two_pi() * a1.cross(a2) / v;
        break;
      default: break;
    }
  }

  /*std::cout << "a1 = " << a1.transpose() << "\n"; 
  std::cout << "a2 = " << a2.transpose() << "\n"; 
  std::cout << "a3 = " << a3.transpose() << "\n"; 
  std::cout << "b1 = " << b1.transpose() << "\n"; 
  std::cout << "b2 = " << b2.transpose() << "\n"; 
  std::cout << "b3 = " << b3.transpose() << "\n"; 
  getchar();*/

  // antiperiodic boundary condition
  Vector3d antipb_shift(0.0, 0.0, 0.0);
  if (lattice.bc1_periodicity()==bc::antiperiodic) antipb_shift(0) = 0.5/lattice.size1();
  if (lattice.bc2_periodicity()==bc::antiperiodic) antipb_shift(1) = 0.5/lattice.size2();
  if (lattice.bc3_periodicity()==bc::antiperiodic) antipb_shift(2) = 0.5/lattice.size3();

  // k-points & translation vectors
  double x1, x2, x3;
  Vector3i n = {0,0,0};
  Vector3i m = {-lattice.size1()/2, -lattice.size2()/2, -lattice.size3()/2};
  this->clear();
  num_kpoint_ = lattice.num_unitcells();
  for (unsigned i=0; i<num_kpoint_; i++) {
    x1 = static_cast<double>(m(0)+n(0))/lattice.size1() + antipb_shift(0);
    x2 = static_cast<double>(m(1)+n(1))/lattice.size2() + antipb_shift(1);
    x3 = static_cast<double>(m(2)+n(2))/lattice.size3() + antipb_shift(2);
    this->push_back(x1*b1 + x2*b2 + x3*b3);
    //kpoints[i] = x1 * b1 + x2 * b2 + x3 * b3;
    //std::cout << i << ": " << kpoints[i](0) << " " << kpoints[i](1) << " " << kpoints[i](2) << "\n";
    //translation_vectors.push_back(n);
    n = lattice.get_next_bravindex(n);
  }
}

/*
kpoint BlochBasis::mesh_nb_dir1(const unsigned& k) const
{
  return operator[](k) + b1/L1_;
}

kpoint BlochBasis::mesh_nb_dir2(const unsigned& k) const
{
  return operator[](k) + b2/L2_;
}

kpoint BlochBasis::mesh_nb_dir3(const unsigned& k) const
{
  return operator[](k) + b3/L3_;
}
*/

// In the k-mesh, the NN point along dir-1
void BlochBasis::gen_mesh_neighbors(const lattice::Lattice& lattice)
{
  int L1 = lattice.size1();
  int L2 = lattice.size2();
  int L3 = lattice.size3();
  std::vector<std::vector<std::vector<int> > > k_points;
  k_points.resize(L1);
  for (int i=0; i<L1; ++i) k_points[i].resize(L2);
  for (int i=0; i<L1; ++i) {
    for (int j=0; j<L2; ++j) {
      k_points[i][j].resize(L3);
    }
  }
  Vector3i n = {0,0,0};
  for (unsigned k=0; k<num_kpoint_; k++) {
    k_points[n(0)][n(1)][n(2)] = k;
    n = lattice.get_next_bravindex(n);
  }

  nn_list_.resize(num_kpoint_);
  for (int k=0; k<num_kpoint_; k++) nn_list_[k].resize(6);

  int nn;
  for (int m=0; m<L3; ++m) {
    for (int j=0; j<L2; ++j) {
      for (int i=0; i<L1; ++i) {
        int k = k_points[i][j][m];
        // right & left nn
        nn = (i+1 < L1) ? i+1 : 0;
        //std::cout << k << " " << nn <<" "<<m << "------HERE---------\n";
        nn_list_[k][0] = k_points[nn][j][m];
        nn = (i-1 >= 0) ? i-1 : L1-1;
        nn_list_[k][1] = k_points[nn][j][m];
        // top & bottom nn
        nn = (j+1 < L2) ? j+1 : 0;
        nn_list_[k][2] = k_points[i][nn][m];
        nn = (j-1 >= 0) ? j-1 : L2-1;
        nn_list_[k][3] = k_points[i][nn][m];
        // front & back 
        nn = (m+1 < L3) ? m+1 : 0;
        nn_list_[k][4] = k_points[i][j][nn];
        nn = (m-1 >= 0) ? m-1 : L3-1;
        nn_list_[k][5] = k_points[i][j][nn];
      }
    }
  }
  // check
  /*for (int k=0; k<num_kpoint_; k++) {
    std::cout<<"nnx, -nnx = "<< nn_list_[k][0]<<", "<<nn_list_[k][1] << "\n";
    std::cout<<"nny, -nny = "<< nn_list_[k][2]<<", "<<nn_list_[k][3] << "\n";
    std::cout<<"nnz, -nnz = "<< nn_list_[k][4]<<", "<<nn_list_[k][5] << "\n";
    getchar();
  }*/
}

void BlochBasis::make_subspace_basis(const lattice::LatticeGraph& graph)
{
  int orbitals_per_cell = 0;
  unsigned n = graph.lattice().num_basis_sites();
  for (unsigned i=0; i<n; ++i) {
    basis_state s = graph.site(i);
    unsigned uid = graph.site_uid(i);
    if (s != graph.site(uid))
      throw std::logic_error("blochbasis::make_site_basis: unexpected graph property.");
    orbitals_per_cell += graph.site_dim(i);
  }
  subspace_dimension_ = n + orbitals_per_cell;
  subspace_basis_.resize(subspace_dimension_);
  for (unsigned i=0; i<subspace_dimension_; ++i) {
    subspace_basis_[i] = i; 
  }
  null_idx_ = subspace_basis_.size();
  // index of the 'representative state' of a site
  representative_state_idx_.resize(graph.num_sites());
  //translation_vectors_.resize(graph.num_sites());
  for (auto& idx : representative_state_idx_) idx = null_idx_;
  for (unsigned i=0; i<graph.num_sites(); ++i) {
    basis_state s = graph.site(i);
    unsigned uid = graph.site_uid(i);
    representative_state_idx_[s] = uid;
    //translation_vectors_[s] = graph.site_cellcord(i);
    //std::cout << translation_vectors_[s] << "\n"; getchar();
  }

  /*
  subspace_dimension_ = graph.lattice().num_basis_sites();
  subspace_basis_.resize(subspace_dimension_);
  for (unsigned i=0; i<subspace_dimension_; ++i) {
    basis_state s = graph.site(i);
    unsigned uid = graph.site_uid(i);
    if (s != graph.site(uid))
      throw std::logic_error("blochbasis::make_site_basis: unexpected graph property.");
    subspace_basis_[i] = s; 
  }
  null_idx_ = subspace_basis_.size();
  // index of the 'representative state' of a site
  representative_state_idx_.resize(graph.num_sites());
  //translation_vectors_.resize(graph.num_sites());
  for (auto& idx : representative_state_idx_) idx = null_idx_;
  for (unsigned i=0; i<graph.num_sites(); ++i) {
    basis_state s = graph.site(i);
    unsigned uid = graph.site_uid(i);
    representative_state_idx_[s] = uid;
    //translation_vectors_[s] = graph.site_cellcord(i);
    //std::cout << translation_vectors_[s] << "\n"; getchar();
  }
  */
}


#ifdef ON
int BlochBasis::construct(const lattice::Lattice& lattice, const lattice::LatticeGraph& graph)
{
  // reset
  kpoints.resize(lattice.num_unitcells());
  //translation_vectors.clear();
  basis_states.resize(lattice.num_basis_sites());
  state_indices.clear();
  //null_state = lattice.num_sites();
  make_bloch_vectors(lattice);
  make_site_basis(lattice, graph);

  return 0;
}


void BlochBasis::make_site_basis(const lattice::Lattice& lattice, const lattice::LatticeGraph& graph)
{
  /*
  * The basis states are ! combination of site states connected by lattice 
  * translational symmetry. Among these states in a 'bloch cycle', only the state
  * (site) with smallest 'id' (representative state=rs) is stored. 
  * The basis states are to be regarded as: 
  *   |s> = 1/sqrt(N) \sum_R e^{ik.R} |rs> 
  */

  // basis states & indexing
  lattice::vertex_iterator vi, vi_end;
  basis_state s;
  unsigned i;
  for (i=0; i<basis_states.size(); ++i) {
    s = graph.vertex(i);
    //std::cout << graph.vertex_uid(i) << "\n";
    if (graph.vertex_id(s) != i) {
      throw std::logic_error("*error: blochbasis:: unexpected graph property");
    }
    basis_states[i] = s;
    state_indices.insert({s,i});
    //std::cout << state_index[s] << "\n";
  }
  null_index = basis_states.size();
}

basis_state BlochBasis::representative_state(const basis_state& s, const lattice::LatticeGraph& graph,
    Vector3d& R) const
{
  unsigned rs = graph.vertex_uid(s);
  R = graph.vertex_cellcord(s);
  return graph.vertex(rs);
}

#endif

} // end namespace basis
