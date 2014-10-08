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
// implied, of the nil authors.
// -------------------------------------------------------------------

#ifndef __nil_parameter_reader_h
#define __nil_parameter_reader_h

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>

namespace nil
{
  
  /**
   * A class that reads a pre-defined set of coefficients from a
   * parameter file. If no parameter file is present, a dummy
   * parameter file is made.
   *
   * @author Toby D. Young 2014
   */
  class ParameterReader 
    : 
    public dealii::Subscriptor
  {
  public:
    
    /**
     * Constructor. Stores a reference to the ParameterHandler object
     * that is passed to it.
     */
    ParameterReader (dealii::ParameterHandler &paramater_handler,
		     const std::string        &parameter_file);
    
    /**
     * Read in the parameters.
     */
    void read_parameters ();
    
  private:
    
    /**
     * Declare parameters that can be used.
     */
    void declare_parameters ();
    
    /**
     * Parse the parameter file to get the number of material
     * ids. @note This function actually does nothing but return the
     * integer 2, which is stupid.
     */
    unsigned int read_material_id ();
    
    /**
     * Reference to the parameter handler.
     */
    dealii::ParameterHandler &prm_handler;
    
    /**
     * Local copy of the parameter file name.
     */
    const std::string prm_file;

  };  
  
} // namespace nil

#endif // __nil_parameter_reader_h
