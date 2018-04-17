/*---------------------------------------------------------------------------
* task.h: Class representing a task to be run.
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-28 10:57:03
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-21 00:10:06
*----------------------------------------------------------------------------*/
#include <iostream>
#include "inputparams.h"
#include "mpi_comm.h"

#ifndef SCHEDULER_WORKER_H
#define SCHEDULER_WORKER_H

namespace scheduler {

class AbstractWorker
{
public:
  virtual ~AbstractWorker() {};
  virtual int start(const input::Parameters& p) = 0; // start the runs
  virtual int run(const input::Parameters& p) = 0; // run 
  virtual int run(const input::Parameters& parms, const mpi_communicator& mpi_comm) = 0;
  virtual void finish(void) = 0; // mark as finished
  virtual void halt(void) = 0; // halt all runs, simulation is finished        
};

class Worker : public AbstractWorker
{
public:
  Worker();
  virtual ~Worker();
  //virtual void construct(); // needs to be called to finish construction
  int start(const input::Parameters& p) override; // start simulation
  int run(const input::Parameters& p) override; // run 
  int run(const input::Parameters& parms, const mpi_communicator& mpi_comm) override;
  virtual void dostep(void)=0; // do a step
  // bool started() const { return started_;}
  void halt(void) override;
protected:
  bool finished_;
private:
  bool started_; // is the task running?
};


} // end namespace input

#endif
