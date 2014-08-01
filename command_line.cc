
#include "command_line.h"

#include <fstream>
#include <iostream>

namespace nil
{
  
  CommandLine::CommandLine () 
    :
    found_prm_file (false)
  {}
  
  CommandLine::~CommandLine () {}
  
  void CommandLine::parse_command_line (const int    argc,
					char *const *argv)
  {
    // Blindly read in the command line parameters.
    for (int i=0; i<argc; ++i)
      args.push_back (argv[i]);
    
    // Keep reading command line arguements until the number of
    // command line arguements reduces to zero.
    while (args.size ())
      {
	
	// See if there is a parameter file to use.
	if ((!found_prm_file)                       && 
	    ((args.front () == std::string ("-pf")) ||
	     (args.front () == std::string ("--parameter_file"))))
	  {
	    // Get rid of the command...
	    args.pop_front ();
	    
	    // and read in the data.
	    runtime_parameters.prm_file = args.front ();
	    args.pop_front ();
	    found_prm_file = true;
	  }
	
	// Otherwise, give up.
	else
	  {
	    break;
	  }
      }
  }
  
} // namespace nil
