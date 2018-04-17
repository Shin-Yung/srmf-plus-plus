/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-04-17 23:13:06
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-01 22:49:39
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef MPI_COMM_H
#define MPI_COMM_H

#include <iostream>
#include <list>
#include <exception>
#ifdef HAVE_BOOST_MPI
  #include <boost/mpi/environment.hpp>
  #include <boost/mpi/communicator.hpp>
#endif

namespace scheduler {

enum {MP_make_task, MP_task_params, MP_run_task, MP_task_finished,
  MP_quit_tasks};

//const int MP_task_params = 0;
//const int MP_run_task = 1;
//const int MP_quit_tasks = 2;

#ifdef HAVE_BOOST_MPI

using mpi_status = boost::mpi::status;
using plist = std::list<int>;
//  using mpi_environment = boost::mpi::environment;
//  using mpi_communicator = boost::mpi::communicator;
class mpi_environment : public boost::mpi::environment
{
public:
  using environment::environment;
  ~mpi_environment() {}
};

class mpi_communicator : public boost::mpi::communicator
{
public:
  using communicator::communicator;
  mpi_communicator();
  ~mpi_communicator() {}
  bool is_master(void) const { return rank()==0; }
  int master(void) const { return 0; }
  const plist& slave_procs(void) const { return slave_procs_; }
private:
  plist slave_procs_;
};

#else

using plist = std::list<int>;

class mpi_status
{
public:
  mpi_status() {}
  ~mpi_status() {}
  int tag(void) const { throw_exception(); return 0; }
  int source(void) const { throw_exception(); return 0; }
private:
  void throw_exception(void) const
  {
    throw std::logic_error("** mpi_communicator:: not an mpi program");
  }
};

class mpi_environment
{
public:
  mpi_environment() {}
  ~mpi_environment() {}
};

class mpi_communicator
{
public:
  mpi_communicator() {}
  ~mpi_communicator() {}
  int rank(void) const { return 0; }
  int size(void) const { return 0; }
  bool is_master(void) const { return rank()==0; }
  int master(void) const { return 0; }
  void send(int dest, int tag) const
    { throw_exception(); }
  template<typename T> void send(int dest, int tag, const T & value) const
    { throw_exception(); }
  template<typename T> void isend(int dest, int tag, const T & value) const
    { throw_exception(); }
  void isend(int dest, int tag) const { throw_exception(); }
  mpi_status recv(int dest, int tag) const 
    { throw_exception(); return mpi_status(); }
  template<typename T> mpi_status recv(int dest, int tag, const T & value) const
    { throw_exception(); return mpi_status(); }
  mpi_status probe(int dest=0, int tag=0) const 
    { throw_exception(); return mpi_status(); }
  mpi_status iprobe(int dest=0, int tag=0) const 
    { throw_exception(); return mpi_status(); }
  const plist& slave_procs(void) const { return slave_procs_; }
private:
  plist slave_procs_;
  void throw_exception(void) const
  {
    throw std::logic_error("** mpi_communicator:: not an mpi program");
  }
};

#endif



} // end namespace scheduler

#endif
