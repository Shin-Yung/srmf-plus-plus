/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-19 11:24:03
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-05-03 22:35:45
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "rotor.h"
#include <stdexcept>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/math/tools/roots.hpp>

namespace srmf {


/*
Cluster::Cluster(const cluster_t& type, const int& id, const rotor_basis& basis)
  : type_{type}, id_{id}, dim_{basis.dim()}
{
  switch (type_) {
    case cluster_t::SITE: num_sites_ = 1; break;
    case cluster_t::BOND: num_sites_ = 2; break;
    case cluster_t::CELL: num_sites_ = 2; break;
  }
  ham_.resize(dim_,dim_);
}
void Cluster::solve()
{
}
*/


Rotor::Rotor(const input::Parameters& inputs, const model::Hamiltonian& model, 
  const lattice::LatticeGraph& graph, const SR_Params& srparams)
  : rotor_graph_(graph)
{
  std::cout << "----Rotor::Rotor-------\n";
  // Rotor lattice has only one original lattice unit cell 
  num_sites_ = srparams.num_sites();
  num_bonds_ = srparams.num_bonds();
  for (auto vi=rotor_graph_.sites_begin(); vi!=rotor_graph_.sites_end(); ++vi) {
    rotor_graph_.change_site_dimension(vi, 1);
  }

  // rotor parameters
  int theta_min = inputs.set_value("min_rotor_charge", -2);
  int theta_max = inputs.set_value("max_rotor_charge", +2);
  std::string cluster_type = inputs.set_value("rotor_cluster", "SITE");
  boost::to_upper(cluster_type);
  if (cluster_type=="SITE") cluster_type_ = cluster_t::SITE;
  else if (cluster_type=="BOND") cluster_type_ = cluster_t::BOND;
  else if (cluster_type=="CELL") cluster_type_ = cluster_t::CELL;
  else throw std::range_error("Rotor::Rotor: invalid 'rotor_cluster'");

  // clusters & basis states per cluster
  make_clusters(srparams);
  basis_.construct(sites_per_cluster_,theta_min,theta_max);
  dim_ = basis_.dim();

  // Rotor model (signle site approaximation)
  using namespace model;
  std::string name;
  double defval;
  CouplingConstant cc;
  double U = model.get_parameter_value("U");
  U_half_ = 0.5 * U;
  rotor_model_.init(rotor_graph_.lattice());
  rotor_model_.add_parameter(name="U", defval=U);
  rotor_model_.add_parameter(name="mu", defval=0.0);
  rotor_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
  rotor_model_.add_siteterm(name="e^i_theta", cc="1", op::ciup_dag());
  rotor_model_.add_siteterm(name="hubbard", cc="U", op::hubbard_int());
  // finalize model
  rotor_model_.finalize(rotor_graph_.lattice());

  // storages
  bond_tchi_.resize(num_bonds_);
  site_density_.resize(num_sites_);
  site_mu_.resize(num_sites_);
  site_phi_.resize(num_sites_);
  rotor_mat_.resize(dim_,dim_);
  rotor_mat_.setZero();
  cluster_ham_.resize(dim_, dim_);
  cluster_ham_.setZero();
  cluster_groundstate_.resize(dim_, clusters_.size());
  for (int i=0; i<num_sites_; ++i) {
    site_mu_[i] = 0.0;
    site_phi_[i] = 1.0;
  }

  construct_matrix();
  solve_clusters();
  find_site_density();
  find_site_phi();
}

void Rotor::solve(SR_Params& srparams) 
{
  // set constrained density
  for (int i=0; i<num_sites_; ++i) site_density_(i) = 1.0-srparams.spinon_density(i);
  constrained_density_ = site_density_.sum()/num_sites_;
  std::cout << "cons_density = " << constrained_density_ << "\n";

  // bond hopping parameters
  bond_tchi_ = srparams.bond_tchi();

  // chemical potential
  double mu = solve_for_mu();
  std::cout << "mu = " << mu << "\n";
}

double Rotor::solve_for_mu(void) 
{
  double guess = 0.0;
  double factor = 2.0;
  const boost::uintmax_t maxit = 20; 
  boost::uintmax_t it = maxit;      
  bool is_rising = true;
  boost::math::tools::eps_tolerance<double> tol(3);
  std::pair<double,double> r = boost::math::tools::bracket_and_solve_root(
    [this](double mu) { 
      for (int i=0; i<num_sites_; ++i) site_mu_[i] = mu;
      solve_clusters();
      find_site_density(); 
      return site_density_.sum()/num_sites_-constrained_density_; 
    }, 
    guess, factor, is_rising, tol, it);
  return r.first + 0.5 * (r.second-r.first);
} 

void Rotor::make_clusters(const SR_Params& srparams)
{
  bonds_ = srparams.bonds();
  site_links_ = srparams.site_links();

  clusters_.clear();
  switch (cluster_type_) {
    case cluster_t::SITE:
      sites_per_cluster_ = 1;
      clusters_.resize(num_sites_);
      for (unsigned i=0; i<num_sites_; ++i) clusters_[i].push_back(i);
      break;

    case cluster_t::BOND:
      sites_per_cluster_ = 2;
      clusters_.resize(bonds_.size());
      for (int i=0; i<bonds_.size(); ++i) clusters_[i].push_back(i);
      break;

    case cluster_t::CELL:
      sites_per_cluster_ = num_sites_;
      clusters_.resize(1);
      for (int i=0; i<bonds_.size(); ++i) clusters_[0].push_back(i);
      break;
  }
}

void Rotor::solve_clusters(void)
{
  switch (cluster_type_) {
    case cluster_t::SITE:
      for (int i=0; i<clusters_.size(); ++i) {
        cluster_ham_.setZero();
        // diagonal elements
        double mu = site_mu_[i];
        for (const auto& elem: diagonal_elems_) {
          int n = elem.row();
          double ntheta = elem.value();
          cluster_ham_(n,n) = (U_half_ * ntheta - mu) * ntheta;
        }
        // site operators
        auto phi = site_phi_[i];
        for (const auto& elem: siteop_elems_) {
          int m = elem.row();
          int n = elem.col();
          cluster_ham_(m,n) = phi;
          cluster_ham_(n,m) = std::conj(phi);
        }
        // solve the hamiltonian
        eigen_solver_.compute(cluster_ham_);
        // store ground state
        //std::cout << "gndstate = " << eigen_solver_.eigenvectors().col(0) << "\n";
        cluster_groundstate_.col(i) = eigen_solver_.eigenvectors().col(0);
        //std::cout << "cluster = " << i << "\n";
        //std::cout << cluster_ham_ << "\n";
        //getchar();
      }
      break;

    case cluster_t::BOND:
      break;

    case cluster_t::CELL:
      break;
  }
  //std::cout << "groundstate = " << cluster_groundstate_ << "\n";
}

void Rotor::find_site_density(void)
{
  double norm, ntheta, sum;
  switch (cluster_type_) {
    case cluster_t::SITE:
      for (int i=0; i<clusters_.size(); ++i) {
        sum = 0.0;
        for (auto n=0; n<dim_; ++n) {
          norm = std::norm(cluster_groundstate_(n,i));
          ntheta = basis_.apply_ni(n);
          sum += norm * ntheta;
        }
        site_density_[i] = sum;
        std::cout << "site density = " << i << ": " << sum << "\n";
      }
      break;
    case cluster_t::BOND:
      break;
    case cluster_t::CELL:
      break;
  }
}

void Rotor::find_site_phi(void)
{
  rotor_basis::idx_t j; // i
  std::complex<double> sum;
  switch (cluster_type_) {
    case cluster_t::SITE:
      for (int n=0; n<clusters_.size(); ++n) {
        sum = 0.0;
        for (j=0; j<dim_-1; ++j) {
          auto c_j = cluster_groundstate_(j,n);
          /*i = basis_.apply_cidag(j);
          if (i != basis_.null_idx()) {
            auto c_i = cluster_groundstate_(i,n);
            sum += std::conj(c_i)*c_j;
          }*/
          auto c_i = cluster_groundstate_(j+1,n);
          sum += std::conj(c_i)*c_j;
        }
        site_phi_[n] = sum;
        std::cout << "site phi = " << n << ": " << sum << "\n";
      }
      break;
    case cluster_t::BOND:
      break;
    case cluster_t::CELL:
      break;
  }
}

void Rotor::construct_matrix(void) 
{
  // non-zero matrix elements of various operators in the Hamiltonian
  rotor_basis::idx_t i, j;
  int mat_elem;
  unsigned site = 0;

  diagonal_elems_.clear();
  for (i=0; i<dim_; ++i) {
    std::tie(mat_elem, j) = basis_.apply_ni(i, site);
    diagonal_elems_.push_back({i,i,mat_elem});
  }

  siteop_elems_.clear();
  for (auto sterm=rotor_model_.siteterms_begin(); sterm!=rotor_model_.siteterms_end(); ++sterm) {
    if (sterm->qn_operator().id()==model::op_id::ciup_dag) {
      //std::cout << "op = " << sterm->qn_operator().name() << "\n";
      for (i=0; i<dim_; ++i) {
        std::tie(mat_elem, j) = basis_.apply_cidag(i, site);
        if (j != basis_.null_idx()) {
          siteop_elems_.push_back({i,j,1});
          //std::cout << "(i,j)=" << i << " " << j << "\n";
        }
      }
    }
  }
  //for (auto& elem : siteop_elems_) std::cout << elem.row() << " " << elem.col() << "\n";
  bondop_elems_.clear();

  // hamiltonian matrix
  site = 0;
  double mu = 0.0;
  for (const auto& elem : diagonal_elems_) {
    i = elem.row();
    double ntheta = elem.value();
    rotor_mat_(i,i) = (U_half_ * ntheta - mu) * ntheta;
  }

  double phi = 1.0;
  for (const auto& elem : siteop_elems_) {
    rotor_mat_(elem.row(), elem.col()) = -phi;
  }
  //std::cout << "rotor_mat = \n" << rotor_mat_ << "\n";
}


} // end namespace srmf
