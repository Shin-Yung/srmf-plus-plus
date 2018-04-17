/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-04-17 23:49:24
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-01 22:49:42
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./mpi_comm.h"

namespace scheduler {

#ifdef HAVE_BOOST_MPI

mpi_communicator::mpi_communicator(void)
  : boost::mpi::communicator()
{
  for (int i=0; i<this->size(); ++i) {
    if (i != this->rank()) slave_procs_.push_back(i);
  }
}

#endif


} // end namespace input
