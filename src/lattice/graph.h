/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: amedhi
* Date:   2016-03-01 00:11:01
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-03 11:47:02
*----------------------------------------------------------------------------*/
#ifndef GRAPH_H
#define GRAPH_H

#include "lattice.h"
#include <stdexcept>
#include <boost/graph/adjacency_list.hpp>
//#include <Eigen/Dense>

namespace lattice {
//namespace graph {

// Vertex properties
struct VertexProperties {
  unsigned uid; // unitcell id
  unsigned dim; // = num_orbitals
  unsigned type; 
  unsigned stype; // symmetry type
  unsigned atomid; 
  Vector3i bravindex;
  Vector3d coord; 
  Vector3d cell_coord;
};

// Edge properties
struct EdgeProperties {
  unsigned type; 
  unsigned stype; 
  int sign; 
  unsigned vector_id; 
  Vector3d vector; 
};

using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, 
              VertexProperties, EdgeProperties>;


class LatticeGraph : public Graph 
{
public:
  //using vertex_descriptor = boost::graph_traits<Graph>::vertex_descriptor;
  //using edge_descriptor = boost::graph_traits<Graph>::edge_descriptor;
  //using vertex_iterator = boost::graph_traits<Graph>::vertex_iterator;
  //using edge_iterator = boost::graph_traits<Graph>::edge_iterator;
  //using out_edge_iterator = boost::graph_traits<Graph>::out_edge_iterator;
  using site_descriptor = boost::graph_traits<Graph>::vertex_descriptor;
  using site_iterator = boost::graph_traits<Graph>::vertex_iterator;
  using bond_iterator = boost::graph_traits<Graph>::edge_iterator;
  using out_bond_iterator = boost::graph_traits<Graph>::out_edge_iterator;
  using in_bond_iterator = boost::graph_traits<Graph>::in_edge_iterator;
  //using adjacency_iterator = boost::graph_traits<Graph>::adjacency_iterator;

  // ctors
  LatticeGraph() { clear(); }
  LatticeGraph(const input::Parameters& parms); 
  //LatticeGraph(const Lattice& lattice);
  ~LatticeGraph() {}

  // setter functions
  void construct(const input::Parameters& parms); 
  void change_type_value(const site_descriptor& v, const unsigned& new_val) 
    { this->operator[](v).type = new_val; }
  void change_type_value(const site_iterator& vi, const unsigned& new_val) 
    { this->operator[](*vi).type = new_val; }
  void change_site_dimension(const site_descriptor& v, const unsigned& new_val) 
    { this->operator[](v).dim = new_val; }
  void change_site_dimension(const site_iterator& vi, const unsigned& new_val) 
    { this->operator[](*vi).dim = new_val; }
  void change_type_value(const bond_iterator& ei, const unsigned& new_val) 
    { this->operator[](*ei).type = new_val; }
  void change_type_value(const out_bond_iterator& ei, const unsigned& new_val) 
    { this->operator[](*ei).type = new_val; }
  void change_type_value(const in_bond_iterator& ei, const unsigned& new_val) 
    { this->operator[](*ei).type = new_val; }

  // getter functions
  const Lattice& lattice(void) const { return lattice_; }

  // vertex accessors
  const unsigned& num_sites(void) const { return num_vertices_; }
  std::pair<const site_iterator, const site_iterator> sites(void) const 
    { return std::make_pair(vi_begin_, vi_end_); }
  const site_iterator& sites_begin(void) const { return vi_begin_; }
  const site_iterator& sites_end(void) const { return vi_end_; }
  unsigned site(const site_iterator& vi) const { return static_cast<unsigned>(*vi); }
  unsigned site(const site_descriptor& v) const { return static_cast<unsigned>(v); }
  const std::set<unsigned>& site_types(void) const { return vertex_types_set_; } 
  unsigned num_site_types(void) const { return vertex_types_set_.size(); } 
  const unsigned& site_type(const site_iterator& vi) const { return this->operator[](*vi).type; }
  const unsigned& site_type(const site_descriptor& v) const { return this->operator[](v).type; }
  const unsigned& site_dim(const site_iterator& vi) const { return this->operator[](*vi).dim; }
  const unsigned& site_dim(const site_descriptor& v) const { return this->operator[](v).dim; }
  const unsigned& site_uid(const site_iterator& vi) const { return this->operator[](*vi).uid; }
  const unsigned& site_uid(const site_descriptor& v) const { return this->operator[](v).uid; }
  const unsigned& site_uid(const unsigned& i) const { return this->operator[](i).uid; }
  const Vector3d& site_cellcord(const site_iterator& vi) const { return this->operator[](*vi).cell_coord; }
  const Vector3d& site_cellcord(const site_descriptor& v) const { return this->operator[](v).cell_coord; }
  const Vector3i& site_bravindex(const site_iterator& vi) const { return this->operator[](*vi).bravindex; }
  const Vector3i& site_bravindex(const site_descriptor& v) const { return this->operator[](v).bravindex; }
  const site_descriptor translated_site(const site_descriptor& vi, const Vector3i& bravidx) const;

  /*
  vertex_descriptor vertex(const unsigned& i) const { return boost::vertex(i, *this); }
  unsigned vertex_id(const vertex_descriptor& v) const { return vertex_index_map[v]; }
  unsigned vertex_id(const vertex_iterator& vi) const { return vertex_index_map[*vi]; }
  unsigned vertex_uid(const vertex_descriptor& v) const { return this->operator[](v).uid; }
  Vector3d vertex_cellcord(const vertex_descriptor& v) const { return this->operator[](v).cell_coord; }
  //unsigned vertex_uid(const unsigned& i) const { return g[boost::vertex(i, g)].uid; }
  unsigned vertex_type(const vertex_descriptor& v) const { return this->operator[](v).type; }
  unsigned vertex_type(const vertex_iterator& vi) const { return this->operator[](*vi).type; }
  */

  // edge accessors
  const bond_iterator& bonds_begin(void) const { return ei_begin_; }
  const bond_iterator& bonds_end(void) const { return ei_end_; }
  const unsigned& num_bonds(void) const { return num_edges_; }
  unsigned num_bond_types(void) const { return edge_types_set_.size(); } 
  const unsigned& bond_type(const bond_iterator& ei) const { return this->operator[](*ei).type; }
  const unsigned& bond_type(const out_bond_iterator& ei) const { return this->operator[](*ei).type; }
  const unsigned& bond_type(const in_bond_iterator& ei) const { return this->operator[](*ei).type; }
  const int& bond_sign(const bond_iterator& ei) const { return this->operator[](*ei).sign; }
  const int& bond_sign(const in_bond_iterator& ei) const { return this->operator[](*ei).sign; }
  const int& bond_sign(const out_bond_iterator& ei) const { return this->operator[](*ei).sign; }
  const unsigned& vector_id(const bond_iterator& ei) const { return this->operator[](*ei).vector_id; }
  const unsigned& vector_id(const in_bond_iterator& ei) const { return this->operator[](*ei).vector_id; }
  const unsigned& vector_id(const out_bond_iterator& ei) const { return this->operator[](*ei).vector_id; }
  const Vector3d& vector(const bond_iterator& ei) const { return this->operator[](*ei).vector; }
  const Vector3d& vector(const in_bond_iterator& ei) const { return this->operator[](*ei).vector; }
  const Vector3d& vector(const out_bond_iterator& ei) const { return this->operator[](*ei).vector; }
  site_descriptor source(const bond_iterator& ei) const { return boost::source(*ei, *this); }
  site_descriptor source(const out_bond_iterator& ei) const { return boost::source(*ei, *this); }
  site_descriptor source(const in_bond_iterator& ei) const { return boost::source(*ei, *this); }
  site_descriptor target(const bond_iterator& ei) const { return boost::target(*ei, *this); }
  site_descriptor target(const out_bond_iterator& ei) const { return boost::target(*ei, *this); }
  site_descriptor target(const in_bond_iterator& ei) const { return boost::target(*ei, *this); }
  //using out_edges = boost::out_edges;
  std::pair<out_edge_iterator, out_edge_iterator> out_bonds(const site_iterator& vi) const 
    { return boost::out_edges(*vi, *this); }
  std::pair<out_bond_iterator, out_bond_iterator> out_bonds(const site_descriptor& v) const 
    { return boost::out_edges(v, *this); }
  std::pair<in_bond_iterator, in_bond_iterator> in_bonds(const site_iterator& vi) const 
    { return boost::in_edges(*vi, *this); }
  std::pair<in_bond_iterator, in_bond_iterator> in_bonds(const site_descriptor& v) const 
    { return boost::in_edges(v, *this); }
  // friends
private:
  lattice::Lattice lattice_;
  unsigned num_vertices_{0};
  unsigned num_edges_{0};
  vertex_iterator vi_begin_, vi_end_;
  edge_iterator ei_begin_, ei_end_;
  boost::property_map<Graph, boost::vertex_index_t>::type vertex_index_map;
  //boost::property_map<graph::Graph, unsigned VertexProperties::*>::type vertex_uid_map;
  //boost::property_map<graph::Graph, Vector3d VertexProperties::*>::type vertex_cellcord_map;
  std::set<unsigned> vertex_types_set_;
  std::set<unsigned> edge_types_set_;
  void construct_graph(void);
};


//} // end namespace graph
} // end namespace lattice

#endif
