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

#include "../include/nil/parameter_reader.h"

namespace nil
{

  ParameterReader::ParameterReader (dealii::ParameterHandler &paramater_handler,
				    const std::string        &parameter_file)
    :
    prm_handler (paramater_handler),
    prm_file    (parameter_file)
  {}
  
  
  unsigned int ParameterReader::read_material_id ()
  {
    // @todo This implemetation is stupid! What this function really
    // should do is parse the parameter file for the number of
    // material ids and use that to declare the correct number of
    // parameter subsections (one for each id).
    return 2;
  }
  
  
  void ParameterReader::read_parameters ()
  {
    declare_parameters ();
    prm_handler.read_input (prm_file);
  }
  
  
  void ParameterReader::declare_parameters ()
  {
    for (unsigned int i=0; i<read_material_id (); ++i)
      {
	std::string subsection 
	  = "Material id " + dealii::Utilities::int_to_string (i);
	
	prm_handler.enter_subsection (subsection);
	{
	  prm_handler.declare_entry ("First-order elastic coefficients",
				    "-0.230",
				    dealii::Patterns::List (dealii::Patterns::Double (), 1),
				    "A list of the first-order elastic coefficients. " 
				    "Default is zinc-blende GaAs. ");
	  
	  prm_handler.declare_entry ("First-order dielectric coefficients",
				    "-0.230",
				    dealii::Patterns::List (dealii::Patterns::Double (), 1),
				    "A list of the first-order dielectric coefficients. " 
				    "Default is zinc-blende GaAs. ");
	  
	  prm_handler.declare_entry ("First-order piezoelectric coefficients",
				    "-0.230",
				    dealii::Patterns::List (dealii::Patterns::Double (), 1),
				    "A list of the first-order piezoelectric coefficients. " 
				    "Default is zinc-blende GaAs. "
				    "See: PRL 96, 187602 (2006). ");
	  
	  prm_handler.declare_entry ("Second-order piezoelectric coefficients",
				    "-0.439, -3.765, -0.492",
				    dealii::Patterns::List (dealii::Patterns::Double (), 1),
				    "A list of the second-order piezoelectric coefficients. "
				    "Default is zinc-blende GaAs. "
				    "See: PRL 96, 187602 (2006). ");
	  
	  prm_handler.declare_entry ("First-order spontaneous polarization coefficients",
				    "0.",
				    dealii::Patterns::List (dealii::Patterns::Double (), 1),
				    "A list of the first-order spontaneous polarization coefficients. " 
				    "Default is zinc-blende GaAs. ");
	  
	  prm_handler.declare_entry ("Bravais lattice coefficients",
				    "1.",
				    dealii::Patterns::List (dealii::Patterns::Double (), 1),
				    "A list of the Bravais lattice coefficients. ");
	}
	prm_handler.leave_subsection ();
      }
  }


  // void ParameterReader::read_parameters ()
  // {}


} // namespace nil


