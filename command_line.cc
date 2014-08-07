// -------------------------------------------------------------------
// @author Toby D. Young
//
// Copyright 2014 nil authors. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE NAMEPSACE EWALENA AUTHORS ``AS
// IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// NAMESPACE EWALENA AUTHORS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
// OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The views and conclusions contained in the software and
// documentation are those of the authors and should not be
// interpreted as representing official policies, either expressed or
// implied, of the namespace ewalena authors.
// -------------------------------------------------------------------

#include "command_line.h"

#include <fstream>
#include <iostream>

namespace nil
{
  
  CommandLine::CommandLine () 
    :
    found_help     (false),
    found_prm_file (false)
  {}

  
  CommandLine::~CommandLine () 
  {}


  std::string CommandLine::get_prm_file ()
  {
    return runtime_parameters.prm_file;
  }
  

  void CommandLine::parse_command_line (const int    argc,
					char *const *argv)
  {

    // Blindly read in the command line parameters.
    for (int i=1; i<argc; ++i)
      args.push_back (argv[i]);
    
    // Keep reading command line arguements until the number of
    // command line arguements reduces to zero.
    while (args.size ())
      {

	// See if help is needed.
	if ((!found_help) && 
	    (args.front () == std::string ("--help")))
	  {
	    // Get rid of the command...
	    args.pop_front ();
	    
	    // and read in the data.
	    runtime_parameters.prm_file = args.front ();
	    args.pop_front ();
	    found_help = true;
	  }

	// See if there is a parameter file to use.
	if ((!found_prm_file)                       && 
	    ((args.front () == std::string ("-pf")) ||
	     (args.front () == std::string ("--parameter-file"))))
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

    if ((!found_prm_file) 
	||
	(found_help))
      {
	// write a usage message to terminal.
	std::cout << std::endl << std::endl
		  << "Usage: step-0 [option]... [file]..."
		  << std::endl << std::endl
		  << "       --help           write this message and exit. "
		  << std::endl
		  << "  -pf, --parameter-file name of the parameter file. "
		  << std::endl << std::endl;
	
	// and exit nicely.
	exit (0);
      }

  }
  
} // namespace nil
