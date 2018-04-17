/*---------------------------------------------------------------------------
* task.cc: Implementation of the task classes.  
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-28 11:01:37
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-21 00:17:41
*----------------------------------------------------------------------------*/
#include "worker.h"

namespace scheduler {

Worker::Worker() : finished_(false), started_(false)
{
}

Worker::~Worker()
{
}

int Worker::start(const input::Parameters& p)
{
  started_ = true;
  return 0;
}

int Worker::run(const input::Parameters& p)
{
  //if(started() && !finished_)
  //dostep();
  return 0;
}

int Worker::run(const input::Parameters& parms, const mpi_communicator& mpi_comm) 
{
  return 0;
}

void Worker::halt()
{
  started_=false;
}

/*
void Task::finish()
{
  finished_=true;
}

// halt all active runs
*/


} // end namespace input
