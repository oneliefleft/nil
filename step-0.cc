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


// deal.II headers
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/tensor.h>

// Library-based headers.
#include "group_symmetry.h"
#include "piezoelectric_tensor.h"
#include "strain_tensor.h"

// Application-base headers.
#include "command_line.h"

#include <fstream>
#include <iostream>


template <enum nil::GroupSymmetry GroupSymm, typename ValueType = double>
class Step0
{
public:

  // First, are the usual class constructors and destructors.
  Step0 ();
  ~Step0 ();
  
  void run ();
  
  // Coefficients are held in a parameter file, so it will be
  // necessary to read in a file from the command line.
  nil::CommandLine command_line;
  
private:
  
  // First is a list of functions that belong to this class (they are
  // explained later on).
  void get_parameters ();
  void setup_system ();
  void assemble_strain_tensor ();
  void solve ();
  void output_results () const;


  // Record of the number of space dimensions we are using.
  const int dim = 3;

  // Following that we have a list of the tensors that will be used in
  // this calculation. They are, first- and second-order piezoelectric
  // tensors
  nil::PiezoelectricTensor<GroupSymm, 1, ValueType> first_order_piezoelectric_tensor;
  nil::PiezoelectricTensor<GroupSymm, 2, ValueType> second_order_piezoelectric_tensor;

  // and a Green strain tensor.
  nil::StrainTensor<GroupSymm, 1, ValueType> strain_tensor;
  
  // Additionally, lists of coefficients are needed for those tensors
  // that are tensors of empirical moduli
  std::vector<ValueType> first_order_piezoelectric_coefficients;
  std::vector<ValueType> second_order_piezoelectric_coefficients;

  // as well as a list of the size of the Bravais lattice in the
  // relaxed configuration
  std::vector<ValueType> bravais_lattice;

  // and in the `actual' configuration
  std::vector<ValueType> actual_bravais_lattice;

  // Specific to this calculation, we need a list of the stretches
  // that will be applied to the primitive cell
  std::vector<ValueType> stretch_tensor_range;

  // and the step-size (or increment)
  double increment;

  // Then we need an object to hold various run-time parameters that
  // are specified in an "prm file".
  dealii::ParameterHandler parameters;
};


// The constructor is typically borning...
template <enum nil::GroupSymmetry GroupSymm, typename ValueType>
Step0<GroupSymm, ValueType>::Step0 ()
{}


// as is the destructor.
template <enum nil::GroupSymmetry GroupSymm, typename ValueType>
Step0<GroupSymm, ValueType>::~Step0 ()
{}


// The first step is to obtain a complete list of the coefficients
// needed for this calculation. First comes a declaration of the
// entries expected to be find in the parameter file and then they are
// read into the object parameters.
template <enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
Step0<GroupSymm, ValueType>::get_parameters ()
{
  // First declare the parameters that are expected to be found.
  parameters.declare_entry ("First-order piezoelectric coefficients",
			    "-0.230",
			    dealii::Patterns::List (dealii::Patterns::Double (), 1),
			    "A list of the first-order piezoelectric coefficients. " 
			    "Default is zinc-blende GaAs. "
			    "See: PRL 96, 187602 (2006). ");
  
  parameters.declare_entry ("Second-order piezoelectric coefficients",
			    "-0.439, -3.765, -0.492",
			    dealii::Patterns::List (dealii::Patterns::Double (), 1),
			    "A list of the second-order piezoelectric coefficients. "
			    "Default is zinc-blende GaAs. "
			    "See: PRL 96, 187602 (2006). ");

  parameters.declare_entry ("Bravais lattice dimensions",
			    "1., 1., 1.",
			    dealii::Patterns::List (dealii::Patterns::Double (), 1),
			    "A list of the Bravais lattice dimensions. ");

  parameters.declare_entry ("Stretch tensor range",
			    "-0.2, 0.2, 0., 0., 0., 0.",
			    dealii::Patterns::List (dealii::Patterns::Double (), 1),
			    "A list of diagonal stretch tensor increments. "
			    "These come in the order -x, +x, -y, +y, -z, +z");

  parameters.declare_entry ("Stretch tensor increment",
			    "0.05",
			    dealii::Patterns::Double (),
			    "The size of diagonal stretch tensor increment.");

  parameters.read_input (command_line.get_prm_file ());

  first_order_piezoelectric_coefficients = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("First-order piezoelectric coefficients")));

  second_order_piezoelectric_coefficients = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("Second-order piezoelectric coefficients")));

  bravais_lattice = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("Bravais lattice dimensions")));

  stretch_tensor_range = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("Stretch tensor range")));

  increment = 
    dealii::Utilities::string_to_double (parameters.get ("Stretch tensor increment"));

  // Some checks
  AssertThrow (bravais_lattice.size ()==3,
	       dealii::ExcMessage ("The Bravais lattice is required to have three components."));

  AssertThrow (stretch_tensor_range.size ()==6,
	       dealii::ExcMessage ("The stretch tensor range is required to have six components."));

  AssertThrow (increment>=0,
	       dealii::ExcMessage ("The stretch tensor increment is required to be >= zero."));
}


// Next initialise all of the objects we are going to use. 
template <enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
Step0<GroupSymm, ValueType>::setup_system ()
{
  // Initialise the first- and second-order piezoelectric tensors...
  first_order_piezoelectric_tensor.reinit ();
  second_order_piezoelectric_tensor.reinit ();

  // and distribute the coefficients.
  first_order_piezoelectric_tensor.distribute_coefficients (first_order_piezoelectric_coefficients);
  second_order_piezoelectric_tensor.distribute_coefficients (second_order_piezoelectric_coefficients);

  // Make sure the `actual' Bravais lattice has the right number of
  // components (3)
  actual_bravais_lattice.resize (3);

  // and fill in values for the initial actual bravais lattice
  // size. @note Because of the way the subroutine assemble_system
  // works, the initial actual Bravais lattice size is equalt to the
  // equilibrium lattice size, plus the lower range (often negative)
  // minus the increment.
  for (unsigned int i=0; i<dim; ++i)
    actual_bravais_lattice[i] = 
      bravais_lattice[i]+stretch_tensor_range[2*i];
}


// Assemble the strain tensor for this system. @note This currently
// works for Green's strain only.
template <enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
Step0<GroupSymm, ValueType>::assemble_strain_tensor ()
{
  // Reinitialise the green strain tensor 
  strain_tensor.reinit ();

  // and fill in the correct values for the diagonal part @note this
  // runs over x-axis only!
  for (unsigned int i=0; i<1; ++i)
    strain_tensor[i][i] = 
      0.5 * (pow ((bravais_lattice[i]/actual_bravais_lattice[i]), 2) - 1);
}


// This routine calculates the quantities of interest for a given
// system.
template <enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
Step0<GroupSymm, ValueType>::solve ()
{
  
}


// Simply output the results to file(s).
template <enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
Step0<GroupSymm, ValueType>::output_results () const
{
  
}


// This is the run function, which wraps all of the above into a
// single logical routine.
template <enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
Step0<GroupSymm, ValueType>::run ()
{
  // First find the parameters need for this calculation
  get_parameters ();

  // and output some data, on the piezoelectric part
  std::cout << "Piezoelectric tensor order:    "
	    << first_order_piezoelectric_tensor.order ()
	    << std::endl
	    << "   Number of space dimensions: "
	    << first_order_piezoelectric_tensor.dim ()
	    << std::endl
	    << "   Rank:                       "
	    << first_order_piezoelectric_tensor.rank ()
	    << std::endl
	    << "   Group symmetry:             "
	    << first_order_piezoelectric_tensor.group_symmetry ()
	    << std::endl
	    << "   Number of coefficients:     "
	    << first_order_piezoelectric_coefficients.size ()
	    << std::endl
	    << "   Values of coefficients:     ";
  for (unsigned int i=0; i<first_order_piezoelectric_coefficients.size (); ++i)
    std::cout << first_order_piezoelectric_coefficients[i] << " ";
  std::cout << std::endl
	    << "Piezoelectric tensor order:    "
	    << second_order_piezoelectric_tensor.order ()
	    << std::endl
	    << "   Number of space dimensions: "
	    << second_order_piezoelectric_tensor.dim ()
	    << std::endl
	    << "   Rank:                       "
	    << second_order_piezoelectric_tensor.rank ()
	    << std::endl
	    << "   Group symmetry:             "
	    << first_order_piezoelectric_tensor.group_symmetry ()
	    << std::endl
	    << "   Number of coefficients:     "
	    << second_order_piezoelectric_coefficients.size ()
	    << std::endl
	    << "   Values of coefficients:     ";
  for (unsigned int i=0; i<second_order_piezoelectric_coefficients.size (); ++i)
    std::cout << second_order_piezoelectric_coefficients[i] << " ";
  std::cout << std::endl;
  
  // and then on the lattice part connected to Green's strain
  std::cout << "Bravais lattice:               "
	    << std::endl
	    << "   Reference dimensions:       ";
  for (unsigned int i=0; i<bravais_lattice.size (); ++i)
    std::cout << bravais_lattice[i] << " ";
  std::cout << std::endl;

  // Having done that now we want to start applying an incremental
  // strain pattern. This is done by taking the ratio of the
  // equilibrium crystal parameters over the `actual' crystal
  // parameters.

  // First, setup the tensors of piezoelectric coefficients
  setup_system ();

  std::cout << std::endl
	    << "----------------------------------------------------"
	    << std::endl;

  // compute the maximum x-range
  const double max_x_range = bravais_lattice[0]+stretch_tensor_range[1]+increment;

  // output some usefull data
  std::cout << "Range x-axis:                  " 
	    << actual_bravais_lattice[0]
	    << " to "
	    << max_x_range-increment
	    << std::endl;

  // Then run over all requested stretches @note As in the above, only
  // run over x-axis for now.
  while (actual_bravais_lattice[0]<max_x_range)
    {
      // assemble the strain tensor for the given increment.
      assemble_strain_tensor ();    
      
      std::cout << "   Actual Bravais lattice:     ";
      for (unsigned int i=0; i<bravais_lattice.size (); ++i)
	std::cout << actual_bravais_lattice[i] << " ";
      std::cout << std::endl
		<< "   Strain tensor:              " << strain_tensor
		<< std::endl;
      
      // Solve for the wanted output
      solve ();
      
      // and push back the results into a table.
      
      
      // Update actual Bravais lattice.
      for (unsigned int i=0; i<dim; ++i)
	actual_bravais_lattice[i] += increment;
    }
  
}


int main (int argc, char **argv)
{

  try
    {
      // Set deal.II's logger depth.
      dealii::deallog.depth_console (0);

      // Initialise the problem.
      Step0<nil::GroupSymmetry::ZincBlende> problem;

      // Parse the command line and input file.
      problem.command_line.parse_command_line (argc, argv); 
      
      // Run the problem.
      problem.run ();
    }

  // ...and if this should fail, try to gather as much information as
  // possible. Specifically, if the exception that was thrown is an
  // object of a class that is derived from the C++ standard class
  // <code>exception</code>.
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }

  // If the exception that was thrown somewhere was not an object of a
  // class derived from the standard <code>exception</code> class,
  // then nothing cane be done at all - simply print an error message
  // and exit.
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  
  // At this point, the application performed as was expected - return
  // without error.
  return 0;
}
