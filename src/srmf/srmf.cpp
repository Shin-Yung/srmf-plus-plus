/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-22 13:00:28
*----------------------------------------------------------------------------*/
#include <iomanip>
#include <fstream>
#include "srmf.h"

namespace srmf {

SRMF::SRMF(const input::Parameters& inputs) 
  : graph_(inputs) 
  , model_(inputs, graph_.lattice())
  , mf_model_(model_,graph_)
  , blochbasis_(graph_)
{
  // outputs
  int info;
  need_chern_number_ = inputs.set_value("chern_number",false,info);
  need_ebands_full_ = inputs.set_value("ebands_full",false,info);
  need_ebands_symm_ = inputs.set_value("ebands_symm",false,info);
  need_band_gap_ = inputs.set_value("band_gap",false,info);

  num_kpoints_ = blochbasis_.num_kpoints();
  kblock_dim_ = blochbasis_.subspace_dimension();
  if (need_chern_number_) {
    blochbasis_.gen_mesh_neighbors(graph_.lattice());
  }

  if (need_ebands_symm_) {
    // High symmetry BZ points for FCC lattice
    Vector3d Gamma = Vector3d(0,0,0);
    Vector3d X = 0.5*(blochbasis_.vector_b2()+blochbasis_.vector_b3());
    Vector3d W = (0.25*blochbasis_.vector_b1()+0.50*blochbasis_.vector_b2()
                     +0.75*blochbasis_.vector_b3());
    Vector3d L = 0.5*(blochbasis_.vector_b1()+blochbasis_.vector_b2()
                     +blochbasis_.vector_b3());
    Vector3d K = (3.0/8*blochbasis_.vector_b1()+3.0/8*blochbasis_.vector_b2()
                     +3.0/4*blochbasis_.vector_b3());
    /*
    std::cout << "W = " << W.transpose() << "\n"; 
    std::cout << "L = " << L.transpose() << "\n"; 
    std::cout << "K = " << K.transpose() << "\n"; 
    getchar();
    */

    int N = 100;
    Vector3d step = (X-Gamma)/N;
    for (int i=0; i<N; ++i) symm_line_.push_back(Gamma+i*step);
    step = (W-X)/N;
    for (int i=0; i<N; ++i) symm_line_.push_back(X+i*step);
    step = (L-W)/N;
    for (int i=0; i<N; ++i) symm_line_.push_back(W+i*step);
    step = (Gamma-L)/N;
    for (int i=0; i<N; ++i) symm_line_.push_back(L+i*step);
    step = (K-Gamma)/N;
    for (int i=0; i<N; ++i) symm_line_.push_back(Gamma+i*step);
    step = (X-K)/N;
    for (int i=0; i<N; ++i) symm_line_.push_back(K+i*step);
  }
}

int SRMF::run(const input::Parameters& inputs) 
{
  //std::cout << " Exiting at SRMF::SRMF\n";
  //std::exit(0);


  // Chern number
  if (need_chern_number_) {
    compute_chern_number();
  }
  // Band gap
  if (need_band_gap_) {
    compute_band_gap();
  }

  // Band structure
  //std::cout << "****SRMF::run:: returning****\n"; return 0;
  //std::ios state(NULL);
  //state.copyfmt(std::cout);
  //std::cout << std::scientific << std::uppercase << std::setprecision(6) << std::right;
  std::ofstream of("res_energy.txt");
  of << std::scientific << std::uppercase << std::setprecision(6) << std::right;
  if (need_ebands_full_) {
    for (unsigned k=0; k<num_kpoints_; ++k) {
      //std::cout << k << " of " << symm_line_.size() << "\n";
      Vector3d kvec = blochbasis_.kvector(k);
      mf_model_.construct_kspace_block(kvec);
      //std::cout << mf_model_.quadratic_spinup_block() << "\n"; std::getchar();
      es_k_up_.compute(mf_model_.quadratic_spinup_block(), Eigen::EigenvaluesOnly);
      //if (k%graph_.lattice().size1()==0) of << "\n";
      of << std::setw(6) << k; 
      of << std::setw(14) << kvec(0) << std::setw(14) << kvec(1); 
      of << std::setw(14) << es_k_up_.eigenvalues().transpose() << "\n";
    }
    //of << "\n\n"; 
    of << "\n"; 
  } 
  if (need_ebands_symm_) {
    for (unsigned k=0; k<symm_line_.size(); ++k) {
      //std::cout << k << " of " << symm_line_.size() << "\n";
      Vector3d kvec = symm_line_[k];
      mf_model_.construct_kspace_block(kvec);
      //std::cout << mf_model_.quadratic_spinup_block() << "\n"; getchar();
      es_k_up_.compute(mf_model_.quadratic_spinup_block(), Eigen::EigenvaluesOnly);
      //ek.insert(ek.end(),es_k_up.eigenvalues().data(),
      //    es_k_up.eigenvalues().data()+kblock_dim_);
      //if (k%graph_.lattice().size1()==0) of << "\n";
      of << std::setw(6) << k; 
      of << std::setw(14) << kvec(0) << std::setw(14) << kvec(1); 
      of << std::setw(14) << es_k_up_.eigenvalues().transpose() << "\n";
    }
    of << "\n\n"; 
  }
  of.close();
  //std::sort(ek.begin(),ek.end());
  //std::cout.copyfmt(state);
  return 0;
}

int SRMF::compute_chern_number(void) 
{
  //using xy_pair = std::pair<std::complex<double>,std::complex<double> >;
  //std::vector<xy_pair> BerryConnection(num_kpoints_);
  //ComplexVector psi_k1(kblock_dim_);
  //ComplexVector psi_k2(kblock_dim_);
  //int band_idx = 0;

  std::vector<ComplexVector> BerryConnection_dir1(num_kpoints_);
  std::vector<ComplexVector> BerryConnection_dir2(num_kpoints_);
  for (int k=0; k<num_kpoints_; ++k) {
    BerryConnection_dir1[k].resize(kblock_dim_);
    BerryConnection_dir2[k].resize(kblock_dim_);
  }
  ComplexMatrix psi_k(kblock_dim_,kblock_dim_);
  for (unsigned k=0; k<num_kpoints_; ++k) {
    Vector3d kv1 = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kv1);
    es_k_up_.compute(mf_model_.quadratic_spinup_block());
    //psi_k1 = es_k_up_.eigenvectors().col(band_idx);
    psi_k = es_k_up_.eigenvectors();

    int k_nn = blochbasis_.mesh_nn_xp(k);
    Vector3d kv2 = blochbasis_.kvector(k_nn);
    mf_model_.construct_kspace_block(kv2);
    es_k_up_.compute(mf_model_.quadratic_spinup_block());
    /*
    psi_k2 = es_k_up_.eigenvectors().col(band_idx);
    auto x = psi_k1.dot(psi_k2);
    BerryConnection[k].first  = x/std::abs(x);
    */
    // Berry Connection at 'k' along dir1 for all bands
    for (int n=0; n<kblock_dim_; ++n) {
      auto x = psi_k.col(n).dot(es_k_up_.eigenvectors().col(n));
      BerryConnection_dir1[k][n] = x/std::abs(x);
    }


    k_nn = blochbasis_.mesh_nn_yp(k);
    kv2 = blochbasis_.kvector(k_nn);
    mf_model_.construct_kspace_block(kv2);
    es_k_up_.compute(mf_model_.quadratic_spinup_block());
    /*
    psi_k2 = es_k_up_.eigenvectors().col(band_idx);
    x = psi_k1.dot(psi_k2);
    BerryConnection[k].second = x/std::abs(x); // along dir-2
    */
    // Berry Connection at 'k' along dir2 for all bands
    for (int n=0; n<kblock_dim_; ++n) {
      auto x = psi_k.col(n).dot(es_k_up_.eigenvectors().col(n));
      BerryConnection_dir2[k][n] = x/std::abs(x);
    }
  }

  std::complex<double> U_01, U_12, U_23, U_30;
  //std::complex<double> csum(0.0);
  ComplexVector BerrySum(kblock_dim_);
  BerrySum.setZero();
  for (unsigned k=0; k<num_kpoints_; ++k) {
    /*
    U_01 = BerryConnection[k].first;
    U_30 = std::conj(BerryConnection[k].second);
    int k_nn = blochbasis_.mesh_nn_xp(k);
    U_12 = BerryConnection[k_nn].second;
    k_nn = blochbasis_.mesh_nn_yp(k);
    U_23 = std::conj(BerryConnection[k_nn].first);
    csum += std::log(U_01 * U_12 * U_23 * U_30);
    */
    int k_nn1 = blochbasis_.mesh_nn_xp(k);
    int k_nn2 = blochbasis_.mesh_nn_yp(k);
    for (int n=0; n<kblock_dim_; ++n) {
      U_01 = BerryConnection_dir1[k][n];
      U_30 = std::conj(BerryConnection_dir2[k][n]);
      U_12 = BerryConnection_dir2[k_nn1][n];
      U_23 = std::conj(BerryConnection_dir1[k_nn2][n]);
      BerrySum[n] += std::log(U_01 * U_12 * U_23 * U_30);
    }
  }
  //std::cout << "csum = " << csum/two_pi() <<"\n";
  //int chern_number = std::nearbyint(std::imag(csum/two_pi()));
  std::vector<int> ChernNumber(kblock_dim_);
  int net_chern_number = 0;
  for (int n=0; n<kblock_dim_; ++n) {
    ChernNumber[n] = std::nearbyint(std::imag(BerrySum[n]/two_pi()));
    net_chern_number += ChernNumber[n];
    std::cout << "Chern number of band-"<<n<<" = "<<ChernNumber[n]<<"\n";
  }
  std::cout << "Net Chern number = " << net_chern_number << "\n";

  return 0;
}

int SRMF::compute_band_gap(void) 
{
  std::ofstream fout("res_bandgap.txt");
  if (kblock_dim_<2) {
    fout << "System does not have multiple bands\n";
    fout.close();
    return 0;
  }

  Eigen::ArrayXd band_lo = RealVector::Constant(kblock_dim_,1.E4);
  Eigen::ArrayXd band_hi = RealVector::Constant(kblock_dim_,-1.E4);
  for (unsigned k=0; k<num_kpoints_; ++k) {
    //std::cout << k << " of " << symm_line_.size() << "\n";
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);
    es_k_up_.compute(mf_model_.quadratic_spinup_block(), Eigen::EigenvaluesOnly);
    band_lo = band_lo.min(es_k_up_.eigenvalues().array());
    band_hi = band_hi.max(es_k_up_.eigenvalues().array());
  }
  //std::cout << "lo = " << band_lo.transpose() << "\n";
  //std::cout << "hi = " << band_hi.transpose() << "\n";
  fout << "Band gaps:" << "\n";
  for (int i=0; i<kblock_dim_-1; ++i) {
    fout<<"Eg("<<i<<","<<i+1<<") = "<<std::setw(14)<<band_lo(i+1)-band_hi(i)<<"\n";
  }
  fout << "\nBand widths:" << "\n";
  for (int i=0; i<kblock_dim_; ++i) {
    fout<<"W("<<i<<") = "<<std::setw(14)<<band_hi(i)-band_lo(i)<<"\n";
  }
  fout << "\nFlatness ratio:" << "\n";
  for (int i=1; i<kblock_dim_-1; ++i) {
    double D1 = band_lo(i)-band_hi(i-1);
    double D2 = band_lo(i+1)-band_hi(i);
    double D = std::min(D1,D2);
    double W = band_hi(i)-band_lo(i);
    fout<<"D/W("<<i<<") = "<<std::setw(14)<<D/W<<"\n";
  }

  fout.close();
  return 0;
}

void SRMF::print_copyright(std::ostream& os)
{
  std::cout << "SRMF\n";
}


} // end namespace diag
