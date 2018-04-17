/*---------------------------------------------------------------------------
* master_scheduler: Scheduler running on master
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-09-28 19:51:04
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-01 22:49:32
*----------------------------------------------------------------------------*/
#include <ctime>
#include <string>
#include <sstream>
#include "scheduler.h"

namespace scheduler {

int start(int argc, const char *argv[], const AbstractTask& theTask)
{
  mpi_environment mpi_env;
  mpi_communicator mpi_comm;
  Scheduler* theScheduler;
  if (mpi_comm.is_master()) {
    theScheduler = new MasterScheduler(argc, argv, mpi_comm, theTask);
    //MasterScheduler master_scheduler(argc, argv, theTask);
    //int res = master_scheduler.run();
  }
  else {
    theScheduler = new Scheduler(mpi_comm, theTask);
  }
  int res = theScheduler->run(mpi_comm);
  return res; 
}

MasterScheduler::MasterScheduler(int argc, const char *argv[], 
  const mpi_communicator& mpi_comm, const AbstractTask& theTask)
  : Scheduler()
  , cmdarg(argc, argv)
{
  if (cmdarg.valid()) {
    if (input.read_inputs(cmdarg)) {
      task_size = input.task_size();
      input.init_task_params();
      input.set_task_params(0);
      // construct 'worker' with the first set of task parameters
      if (!cmdarg.have_option(quiet)) 
        std::cout << " starting..." << std::endl;
      theWorker = theTask.make_worker(input.task_params());
      valid_ = true;
      // send to slave schedulers
      if (mpi_comm.size()>1) {
        for (const auto& p : mpi_comm.slave_procs())
          mpi_comm.isend(p, MP_make_task, input.task_params());
      }
    }
  }
}

int MasterScheduler::run(const mpi_communicator& mpi_comm) 
{
  if (!valid_) return -1;
  auto start_time = std::chrono::steady_clock::now();
  //std::cout << " starting..." << std::endl;
  for (unsigned task=0; task<task_size; ++task) {
    //------------ timing -------------
    auto tp = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(tp);
    std::string ts = std::ctime(&t);
    std::string tm = ts.substr(11,8)+" "+ts.substr(8,2)+"-"+ts.substr(4,3)+"-"+ts.substr(20,4);
    if (!cmdarg.have_option(quiet)) 
      std::cout <<" starting task "<<task+1<<" of "<< task_size<<" at "<<tm<<"\n";
    auto tstart_t = std::chrono::steady_clock::now();
    //------------------------------
    input.set_task_params(task);
    //params << pstore(task_id);
    // send to slave schedulers
    if (mpi_comm.size()>1) {
      for (const auto& p : mpi_comm.slave_procs())
        mpi_comm.isend(p, MP_run_task,input.task_params());
      // run own task
      theWorker->run(input.task_params(),mpi_comm);
      // wave for slave process
      for (const auto& p : mpi_comm.slave_procs())
        mpi_comm.recv(p, MP_task_finished);
    }
    else {
      theWorker->run(input.task_params());
    }

    //-------------------------------------
    auto tend_t = std::chrono::steady_clock::now();
    if (!cmdarg.have_option(quiet)) 
      std::cout << " - finished this task in " << elapsed_time(tstart_t, tend_t) << "\n";
    //-------------------------------------
  }
  // finish
  for (const auto& p : mpi_comm.slave_procs()) {
    mpi_comm.isend(p, MP_quit_tasks);
  }
  theWorker->finish();

  auto end_time = std::chrono::steady_clock::now();
  if (!cmdarg.have_option(quiet)) {
    std::cout << " finished all tasks\n"; 
    std::cout << " time: " << elapsed_time(start_time, end_time) << "\n";
    std::cout << " done!\n";
  }
  return 0;
}

std::string MasterScheduler::elapsed_time(const std::chrono::steady_clock::time_point& start_time,
  const std::chrono::steady_clock::time_point& end_time) const
{
  std::ostringstream os;
  auto diff = end_time - start_time;
  auto secs = std::chrono::duration_cast<std::chrono::seconds>(diff).count();
  auto days = secs/86400;
  secs = secs%86400;
  auto hrs = secs/3600;
  secs = secs%3600;
  auto mins = secs/60;
  secs = secs%60;
  if (days > 0) {
    os << days << "d:" << hrs << "h:" << mins << "m";
    return os.str();
  }
  else if (hrs > 0) {
    os << hrs << "h:" << mins << "m:" << secs << "s";
    return os.str();
  }
  else {
    auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
    mins = milli/60000;
    milli = milli%60000;
    secs = milli/1000;
    milli = milli%1000;
    os << mins << "m:" << secs << "s:" << milli << "ms";
    return os.str();
  }
}

} // end namespace scheduler
