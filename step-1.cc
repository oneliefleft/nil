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
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/table_handler.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

// Library-based headers.
#include "group_symmetry.h"
#include "piezoelectric_tensor.h"

// Application-base headers.
#include "command_line.h"

#include <fstream>
#include <iostream>


template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType = double>
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

  void create_coarse_grid (const unsigned int n_refinement_cycles = 0);

  void setup_system ();

  void output_results () const;

  // Following that we have a list of the tensors that will be used in
  // this calculation. They are, first- and second-order piezoelectric
  // tensors
  nil::PiezoelectricTensor<GroupSymm, 1, ValueType> first_order_piezoelectric_tensor;
  nil::PiezoelectricTensor<GroupSymm, 2, ValueType> second_order_piezoelectric_tensor;
 
  // Additionally, lists of coefficients are needed for those tensors
  // that are tensors of empirical moduli
  std::vector<ValueType> first_order_piezoelectric_coefficients;
  std::vector<ValueType> second_order_piezoelectric_coefficients;

  dealii::ConditionalOStream pcout;

  // A parallel distributed triangulation
  dealii::parallel::distributed::Triangulation<dim> triangulation;

  // Then we need an object to hold various run-time parameters that
  // are specified in an "prm file".
  dealii::ParameterHandler parameters;

  // and finally, a table to handle the results.
  dealii::TableHandler output_table;
};


// The constructor is typically borning...
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
Step0<dim, GroupSymm, ValueType>::Step0 ()
  :
  pcout (std::cout, (dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)),

  triangulation (MPI_COMM_WORLD,
                   typename dealii::Triangulation<dim>::MeshSmoothing
                   (dealii::Triangulation<dim>::smoothing_on_refinement |
		    dealii::Triangulation<dim>::smoothing_on_coarsening))
{}


// as is the destructor.
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
Step0<dim, GroupSymm, ValueType>::~Step0 ()
{}


// The first step is to obtain a complete list of the coefficients
// needed for this calculation. First comes a declaration of the
// entries expected to be find in the parameter file and then they are
// read into the object parameters.
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
Step0<dim, GroupSymm, ValueType>::get_parameters ()
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

  // parameters.declare_entry ("Bravais lattice dimensions",
  // 			    "1., 1., 1.",
  // 			    dealii::Patterns::List (dealii::Patterns::Double (), 1),
  // 			    "A list of the Bravais lattice dimensions. ");

  parameters.read_input (command_line.get_prm_file ());

  first_order_piezoelectric_coefficients = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("First-order piezoelectric coefficients")));

  second_order_piezoelectric_coefficients = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("Second-order piezoelectric coefficients")));

  // bravais_lattice = 
  //   dealii::Utilities::string_to_double
  //   (dealii::Utilities::split_string_list (parameters.get ("Bravais lattice dimensions")));

  // Some checks
  // AssertThrow (bravais_lattice.size ()==3,
  // 	       dealii::ExcMessage ("The Bravais lattice is required to have three components."));

}


template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
Step0<dim, GroupSymm, ValueType>::create_coarse_grid (const unsigned int n_refinement_levels_cycles)
{
  dealii::GridGenerator::hyper_rectangle (triangulation,
					  dealii::Point<dim> (-20,-20,-20),
					  dealii::Point<dim> ( 20, 20, 20),
					  true);

  triangulation.refine_global (n_refinement_cycles);
}

// Next initialise all of the objects we are going to use. 
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
Step0<dim, GroupSymm, ValueType>::setup_system ()
{
  // Initialise the first- and second-order piezoelectric tensors...
  first_order_piezoelectric_tensor.reinit ();
  second_order_piezoelectric_tensor.reinit ();

  // and distribute the coefficients.
  first_order_piezoelectric_tensor.distribute_coefficients (first_order_piezoelectric_coefficients);
  second_order_piezoelectric_tensor.distribute_coefficients (second_order_piezoelectric_coefficients);
}




// This is the run function, which wraps all of the above into a
// single logical routine.
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
Step0<dim, GroupSymm, ValueType>::run ()
{

  // First find the parameters need for this calculation
  get_parameters ();
 
  // Then create the coarse grid
  create_coarse_grid (5);

  pcout << "Number of active cells: "
	<< triangulation.n_global_active_cells ()
	<< " (on "
	<< triangulation.n_levels ()
	<< " levels)"
	<< std::endl;
}


int main (int argc, char **argv)
{

  // Initialise MPI as we allways do.
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 
							       dealii::numbers::invalid_unsigned_int);

  try
    {
      // Set deal.II's logger depth.
      dealii::deallog.depth_console (0);

      // Initialise the problem, 3d wurtzite.
      Step0<3, nil::GroupSymmetry::Wurtzite> problem;

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
