/*---------------------------------------------------------------------------
* Scheduler: A class that handles user jobs.  
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-08-17 13:33:19
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-21 00:12:10
*----------------------------------------------------------------------------*/
// File: scheduler.h 
// Definition of the Scheduler class.

#ifndef SCHEDULER_SCHEDULER_H
#define SCHEDULER_SCHEDULER_H

#include <iostream>
#include <chrono>
#include "./task.h"
#include "./cmdargs.h"
#include "./inputparams.h"
#include "./mpi_comm.h"

namespace scheduler {

int start(int argc, const char *argv[], const AbstractTask& theTask);

class Scheduler 
{
public:
  //Scheduler(): simmaster(0) {};
  //Scheduler(Task& theTask) {}
  Scheduler() {}
  Scheduler(const mpi_communicator& mpi_comm, const AbstractTask& theTask);
  ~Scheduler() {} 
  virtual int run(const mpi_communicator& mpi_comm);

protected:
  //Worker theWorker_;
  AbstractWorker* theWorker;
  input::JobInput input;
  unsigned int task_size;
  bool valid_{false};

private:
  //unsigned simmaster;
  //TaskParams parms;
};

class MasterScheduler : public Scheduler
{
public:
  MasterScheduler(int argc, const char *argv[], const mpi_communicator& mpi_comm, 
    const AbstractTask& theTask);
  //MasterScheduler(int argc, const char *argv[], const Task&);
  MasterScheduler() = delete;
  ~MasterScheduler() {}
  int run(const mpi_communicator& mpi_comm) override;

private:
  CommandArg cmdarg;

  std::string elapsed_time(const std::chrono::steady_clock::time_point& start_time,
    const std::chrono::steady_clock::time_point& end_time) const;
};

} // end namespace input

#endif
