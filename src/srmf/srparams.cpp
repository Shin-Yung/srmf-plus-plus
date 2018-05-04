/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-29 21:46:50
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-05-03 10:27:46
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "srparams.h"

namespace srmf {

SR_Params::SR_Params(const input::Parameters& inputs, const lattice::LatticeGraph& graph)
{
  //std::cout << "----SR_Params::SR_Params----\n";
  num_sites_ = graph.lattice().num_basis_sites();
  // store the bonds in a 'unit cell'
  lattice::LatticeGraph::out_edge_iterator ei, ei_end;
  bonds_.clear();
  site_links_.clear();
  site_links_.resize(num_sites_);
  for (unsigned i=0; i<num_sites_; ++i) {
    for (std::tie(ei, ei_end)=graph.out_bonds(i); ei!=ei_end; ++ei) {
      auto type = graph.bond_type(ei);
      auto s = graph.source(ei);
      auto t = graph.target(ei);
      unsigned m = graph.site_uid(s);
      unsigned n = graph.site_uid(t);
      bonds_.push_back({type,m,n});
      // store id of the bond connected to the site
      int id = bonds_.size()-1;
      site_links_[m].push_back({id,false}); // outgoing from 'm'
      site_links_[n].push_back({id,true}); // incoming to 'n'
      //std::cout<<bonds_.back().type()<<": "<< bonds_.back().src()<<" - "<< bonds_.back().tgt()<<"\n";
    }
  }
  num_bonds_ = bonds_.size();

  // storages
  bond_tchi_.resize(num_bonds_);
  sbond_avg_.resize(num_bonds_);
  rbond_avg_.resize(num_bonds_);
  ssite_avg_.resize(num_sites_);
  rsite_avg_.resize(num_sites_);
  spinon_site_density_.resize(num_sites_);

  for (int i=0; i<num_bonds_; ++i) {
    bond_tchi_[i] = 1.0;
    sbond_avg_[i] = 1.0;
    rbond_avg_[i] = 1.0;
  }
  for (int i=0; i<num_sites_; ++i) {
    ssite_avg_[i] = 0.0;
    rsite_avg_[i] = 0.0;
    spinon_site_density_[i] = 1.0;
  }

  // spinon density

}



} // end namespace srmf
