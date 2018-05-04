/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef DIAG_H
#define DIAG_H

#include <iostream>
#include <Eigen/Eigenvalues>
#include "../scheduler/worker.h"
#include "../lattice/constants.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../lattice/blochbasis.h"
#include "../model/model.h"
#include "srparams.h"
#include "rotor.h"
#include "spinon.h"

namespace srmf {

class SRMF : public scheduler::Worker
{
public:
  SRMF(const input::Parameters& parms); 
  ~SRMF() {}
  int start(const input::Parameters& parms) override { return 0; }
  int run(const input::Parameters& parms) override;
  void finish(void) override {} 
  void dostep(void) override {} 
  void halt(void) override {} 
  static void print_copyright(std::ostream& os);
private:
  lattice::LatticeGraph graph_;
  basis::BlochBasis blochbasis_;
  model::Hamiltonian model_;
  SR_Params sr_parms_;
  Spinon spinon_;
  Rotor rotor_;
  unsigned num_kpoints_{1};
  unsigned kblock_dim_{1};
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_k_up_;
  std::vector<Vector3d> symm_line_;
  // outputs
  /*bool need_chern_number_{false};
  bool need_ebands_full_{false};
  bool need_ebands_symm_{false};
  bool need_band_gap_{false};*/

  int compute_chern_number(void);
  int compute_band_gap(void);
};


} // end namespace srmf

#endif
