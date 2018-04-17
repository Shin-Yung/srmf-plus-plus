/*---------------------------------------------------------------------------
* DMRG Project: DMRG using Matrix Product States 
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-08-17 13:33:19
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-18 00:38:02
*----------------------------------------------------------------------------*/
// File: cmdopts.h 
// Class declarations for handling command line arguments

#ifndef SCHEDULER_CMDOPTS_H
#define SCHEDULER_CMDOPTS_H

#include <iostream>
#include <string>
#include "optionparser.h"

namespace scheduler {

// Options index & type
enum OptionIndex {unknown, info, help, quiet, huge, test};
enum OptionType {disable, enable, other};

// Descriptions of the command line options & arguments
const option::Descriptor usage[] =
{
 {unknown, other, "", "", option::Arg::None, "USAGE: program [options] [file]\n\n" "Options:" },
 {info, other, "", "info", option::Arg::None, "  --info  \tPrint build information and exit." },
 {help, other, "", "help", option::Arg::None, "  --help  \tPrint usage and exit." },
 {quiet, enable, "q", "quiet", option::Arg::None, "  --quiet -q  \tEnable quiet mode." },
 {test, enable, "t", "test", option::Arg::None, "  --test -t  \tEnable test run mode." },
 {unknown, other, "", "",    option::Arg::None, "\nExamples:\n"
                                            "  program --quiet --huge file\n"
                                            "  program -qht file\n"
                                            "  program\n" },
 {0,0,0,0,0,0}  // mandatory, to signal end of options 
};

class CommandArg
{
public:
  CommandArg(int argc, const char *argv[]); // constructor
  CommandArg(); // default constructor
  ~CommandArg(); // destructor
  std::string extract_filename(const std::string& path) const;
  bool process_options(void) const;
  bool have_option(enum OptionIndex idx) const { return options[idx]; }
  bool not_valid(void) const {return !valid_;}
  bool valid(void) const {return valid_;}
  std::string filename(void) const { return inputfile; }

private:
  bool valid_;
  int option_count;
  int nonOption_count;
  std::string progname;
  std::string inputfile;
  option::Option* options;
  option::Option* buffer;
};


} // end namespace

#endif
