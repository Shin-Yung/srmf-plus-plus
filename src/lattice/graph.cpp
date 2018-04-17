/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: amedhi
* Date:   2016-03-01 00:11:01
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-03 11:46:59
*----------------------------------------------------------------------------*/
#include "graph.h"

namespace lattice {
//namespace graph {

//using namespace boost;
//LatticeGraph::LatticeGraph(const Lattice& lattice) 
//{
//  construct(lattice);
//}

LatticeGraph::LatticeGraph(const input::Parameters& parms) : lattice_(parms)
{
  construct_graph();
  //std::cout << "-----------Exiting at LatticeGraph::LatticeGraph----------\n";
  //std::exit(0);
} 

void LatticeGraph::construct(const input::Parameters& parms) 
{
  lattice_.construct(parms);
  construct_graph();
} 

void LatticeGraph::construct_graph(void) 
{
  // all the sites and the bonds
  std::vector<Site> sites;
  std::vector<Bond> bonds;
  Unitcell translated_cell;
  Vector3i bravindex(0,0,0);
  for (unsigned i=0; i<lattice_.num_unitcells(); ++i) {
    translated_cell = lattice_.get_translated_cell(bravindex);
    // collect the sites
    for (unsigned n=0; n<translated_cell.num_sites(); ++n) 
      sites.push_back(translated_cell.site(n));
    // collect the bonds
    for (unsigned n=0; n<translated_cell.num_bonds(); ++n) 
      bonds.push_back(translated_cell.bond(n));
    bravindex = lattice_.get_next_bravindex(bravindex);
  }

  // construct the graph
  clear();
  vertex_types_set_.clear();

  // add vertices 
  for (unsigned i=0; i<sites.size(); ++i) {
    vertex_descriptor v = add_vertex(*this);
    // vertex properties
    this->operator[](v).uid = sites[i].uid();
    this->operator[](v).dim = sites[i].num_orbitals();
    this->operator[](v).type = sites[i].type();
    this->operator[](v).stype = sites[i].type();
    this->operator[](v).atomid = sites[i].atomid();
    this->operator[](v).bravindex = sites[i].bravindex();
    this->operator[](v).coord = sites[i].coord();
    this->operator[](v).cell_coord = sites[i].cell_coord();
    //std::cout << " i " << i << " = " << this->operator[](v).dim << "\n";

    // collect type value in a set
    vertex_types_set_.insert(sites[i].type());
  } 
  num_vertices_ = boost::num_vertices(*this);
  boost::tie(vi_begin_, vi_end_) = boost::vertices(*this);
  vertex_index_map = boost::get(boost::vertex_index, *this);
  // is the vertex_descriptor same as its 'position' or 'id'?
  for (unsigned i=0; i<num_vertices_; ++i) {
    if (i != boost::vertex(i, *this)) {
      throw std::logic_error("LatticeGraph::construct_graph: vertex_descriptor assumption wrong");
    }
  } 

  //vertex_uid_map = get(&VertexProperties::uid, g); 
  //vertex_cellcord_map = get(&VertexProperties::cell_coord, g); 

  // add the edges
  edge_types_set_.clear();
  edge_descriptor e;
  bool flag;
  vertex_descriptor source, target;
  for (unsigned i=0; i<bonds.size(); ++i) {
    Bond b = bonds[i];
    // connect the bond, discard if can't be connected (for bonds across open boundaries) 
    if (!lattice_.connect_bond(b, sites)) continue; 
    // add to graph
    source = boost::vertex(b.src_id(), *this);
    target = boost::vertex(b.tgt_id(), *this);
    boost::tie(e, flag) = add_edge(source, target, *this);
    //std::cout << e << "\n";
    // edge properties
    this->operator[](e).type = b.type();
    this->operator[](e).stype = b.type();
    this->operator[](e).sign = b.sign();
    this->operator[](e).vector_id = b.vector_id();
    //this->operator[](e).vector = sites[b.tgt_id()].coord() - sites[b.src_id()].coord();
    this->operator[](e).vector = b.vector();
    //std::cout << "bond: " << vertex_index_map[source] << " --- " << vertex_index_map[target] << "\n";
    //std::cout << bonds[i].src_offset() << "\n";
    //std::cout << bonds[i].tgt_offset() << "\n\n";
    edge_types_set_.insert(b.type());
  }
  num_edges_ = boost::num_edges(*this);
  boost::tie(ei_begin_, ei_end_) = boost::edges(*this);
  //for (unsigned i=0; i<num_edges_; ++i) {
  //  std::cout << boost::edge(i, *this) << "\n"; 
  //}
  //std::cout << "Vertices = " << num_vertices(g) << std::endl;
  //std::cout << "Edges = " << num_edges(g) << std::endl;
  //std::pair<adjacency_iterator,adjacency_iterator> p = boost::adjacent_vertices(0, *this);
  //for (auto it=p.first; it!=p.second;++it) { std::cout << *it << "\n"; }
  //std::pair<in_edge_iterator, in_edge_iterator> po = boost::in_edges(0, *this); 
  //for (auto it=po.first; it!=po.second;++it) { std::cout << *it << "\n"; }
}

const LatticeGraph::site_descriptor 
LatticeGraph::translated_site(const LatticeGraph::site_descriptor& v, 
  const Vector3i& translation_vec) const
{
  return static_cast<site_descriptor>(lattice_.translation_mapped_site(
    site_uid(v), site_bravindex(v), translation_vec)); 
}

//} // end namespace graph
} // end namespace lattice
