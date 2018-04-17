/*---------------------------------------------------------------------------
* DMRG Project: DMRG using Matrix Product States 
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-08-17 13:33:19
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-08 20:57:43
*----------------------------------------------------------------------------*/
// File: cmdargs.cpp 
// Class definitions for handling command line arguments

//#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include "optionparser.h"
#include "cmdargs.h"

namespace scheduler {

// constructor
CommandArg::CommandArg(int argc, const char *argv[]): progname("program"), inputfile("")
{
  if (argc > 0) {
    //boost::filesystem::path p(argv[0]);
    //progname = p.filename().string();
    progname = extract_filename(argv[0]);
    //std::cout << progname << "\n";
    --argc; ++argv; // skip program name for the rest
  }
  // option stats
  option::Stats stats(usage, argc, argv);
  options = new option::Option[stats.options_max];
  buffer = new option::Option[stats.buffer_max];
  // parse the arguments
  option::Parser parse(usage, argc, argv, options, buffer);
  option_count = parse.optionsCount();
  nonOption_count = parse.nonOptionsCount();
  if (nonOption_count) inputfile = parse.nonOption(0);
  valid_ = process_options();
  if (!valid_) inputfile.clear();
}

// default constructor
CommandArg::CommandArg(): progname("program"), inputfile("input.parm")
{
  options = new option::Option[1];
  buffer = new option::Option[1];
  option_count = 0;
}

// desctructor
CommandArg::~CommandArg()
{
  delete[] options;
  delete[] buffer;
}

// member functions
std::string CommandArg::extract_filename(const std::string& path) const
{
  return path.substr(path.find_last_of("/\\") + 1);
}

bool CommandArg::process_options(void) const
{
  // handle a few of the options
  if (options[unknown]) {
    std::cout << "Unknown option: " << options[unknown].name << "\n\n";
    option::printUsage(std::cout, usage);
    return false;
  }

  if (options[help]) {
    option::printUsage(std::cout, usage);
    return false;
  }

  if (options[info]) {
    std::cout << "Build info:" << std::endl;
    return false;
  }

  //option::Option* opt = options[file];
  //if (options[file]) inputfile = options[file].arg;

  return true;
}


} // end namespace
