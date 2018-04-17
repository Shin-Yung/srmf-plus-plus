/*---------------------------------------------------------------------------
* Scheduler: A class that handles user jobs.  
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-08-17 12:44:04
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-01 22:49:36
*----------------------------------------------------------------------------*/
#include <iostream>
#include "scheduler.h"

namespace scheduler {


Scheduler::Scheduler(const mpi_communicator& mpi_comm, const AbstractTask& theTask)
{
  mpi_comm.recv(mpi_comm.master(),MP_make_task,input.task_params());
  theWorker = theTask.make_worker(input.task_params());
  //std::cout << "I am slave " << mpi_comm.rank() << "\n";
  //std::cout << "id = " << input.task_params().task_id() << "\n";
  //std::cout << "size = " << input.task_params().task_size() << "\n";
}

int Scheduler::run(const mpi_communicator& mpi_comm) 
{
  bool task_exist = false;
  while (true) {
    mpi_status msg = mpi_comm.probe();
    switch (msg.tag()) {
      case MP_quit_tasks:
        mpi_comm.recv(msg.source(),msg.tag());
        theWorker->finish();
        return 0;
      case MP_run_task:
        mpi_comm.recv(msg.source(),msg.tag(),input.task_params());
        theWorker->run(input.task_params(),mpi_comm);
        mpi_comm.send(msg.source(),MP_task_finished);
        task_exist = false;
        break;
      default:
        std::cout << "Scheduler::run: unexpected message\n";
        return 1;
    }
  }
  return 0;
}

} // end namespace scheduler
