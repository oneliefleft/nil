// -------------------------------------------------------------------
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

#include <fstream>
#include <iostream>
#include <list>

namespace nil
{
  
  /**
   * This class holds command line information that is given ar
   * runtume.
   */
  class CommandLine
  {
  public:
    
    /**
     * Constructor.
     */
    CommandLine ();
    
    /**
     * Destructor.
     */
    ~CommandLine ();

    /**
     * Return the name of the parameter file as an <code>std::string</code>.
     */
    std::string get_prm_file ();
    
    /**
     * Parse the command line. @note <code>arg<code> and
     * <code>argv</code> are expected to be the C standard.
     */
    void parse_command_line (const int    argc,
			     char *const *argv);
    
  private:
    
    /**
     * A structure that holds the runtime parameters read in from the
     * command line.
     */
    struct RuntimeParameters
    {
      std::string prm_file;  
    } 
    runtime_parameters;
    
    /**
     * A vector to hold command line arguments.
     */
    std::list<std::string> args;

    /**
     * True if help was found on the command line (otherwise false).
     */    
    bool found_help;

    /**
     * True if the parameter file name was found on the command line
     * (otherwise false). @note The short-hand <code>prm</code> is a
     * deal.II standard.
     */
    bool found_prm_file;
    
  }; // class CommandLine
  
} // namespace nil

