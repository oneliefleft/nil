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
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/sparsity_tools.h>

// Library-based headers.
#include "group_symmetry.h"
#include "dielectric_tensor.h"
#include "elastic_tensor.h"
#include "piezoelectric_tensor.h"

// Application-base headers.
#include "command_line.h"

#include <fstream>
#include <iostream>


/**
 * This class can setup, solve, and output the results of the
 * piezoelectric problem.
 *
 * @author Toby D. Young 2012, 2014
 */ 
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType = double>
class PiezoelectricProblem
{
public:

  /**
   * Constructor.
   */
  PiezoelectricProblem ();

  /**
   * Destructor.
   */
  ~PiezoelectricProblem ();
  
  /**
   * Run the problem (a basic useage example)
   */
  void run ();
  
  // Coefficients are held in a parameter file, so it will be
  // necessary to read in a file from the command line.
  nil::CommandLine command_line;
  
private:
  
  // First is a list of functions that belong to this class (they are
  // explained later on).
  void get_parameters ();

  void make_coarse_grid (const unsigned int n_refinement_cycles = 0);
  void make_boundary_constraints ();

  void setup_system ();
  void setup_coefficients ();
  void assemble_system ();

  void output_results () const;

  /**
   * A local copy of the MPI communicator.
   */
  MPI_Comm mpi_communicator;

  // Following that we have a list of the tensors that will be used in
  // this calculation. They are, first- 
  nil::ElasticTensor<GroupSymm, 1, ValueType>       first_order_elastic_tensor;
  nil::DielectricTensor<GroupSymm, 1, ValueType>    first_order_dielectric_tensor;
  nil::PiezoelectricTensor<GroupSymm, 1, ValueType> first_order_piezoelectric_tensor;
  nil::PiezoelectricTensor<GroupSymm, 1, ValueType> first_order_spontaneous_polarization_tensor;

  // and second-order piezoelectric tensors
  nil::PiezoelectricTensor<GroupSymm, 2, ValueType> second_order_piezoelectric_tensor;
 
  // Additionally, lists of coefficients are needed for those tensors
  std::vector<ValueType> first_order_elastic_coefficients;
  std::vector<ValueType> first_order_dielectric_coefficients;
  std::vector<ValueType> first_order_piezoelectric_coefficients;
  std::vector<ValueType> first_order_spontaneous_polarization_coefficients;

  std::vector<ValueType> second_order_piezoelectric_coefficients;

  /**
   * The size of the Bravais lattice.
   */
  std::vector<double>    bravais_lattice;

  /**
   * A I<code>deal.II</code> hack tha outputs to the first processor
   * only (useful for output in parallel ca;culations).
   */
  dealii::ConditionalOStream pcout;

  /**
   * A parallel distributed triangulation.
   */
  dealii::parallel::distributed::Triangulation<dim> triangulation;

  /**
   * The finite element system.
   */
  const dealii::FESystem<dim>                    fe_q;

  /**
   * Handler of degrees of freedom.
   */
  dealii::DoFHandler<dim>                        dof_handler;

  /**
   * A matrix holding all constraints.
   */
  dealii::ConstraintMatrix                       constraints;

  /**
   * An index set of the degrees of freedom locally owned by this
   * processor.
   */
  dealii::IndexSet                               locally_owned_dofs;

  /**
   * An index set of the degrees of freedom?
   */
  dealii::IndexSet                               locally_relevant_dofs;

  /**
   * The system matrix. This has a block structure but, since we do
   * not need the blocks, us a standard MPI sparse matrix.
   */
  dealii::PETScWrappers::MPI::SparseMatrix       system_matrix;

  /**
   * The system rhs. This has a block structure but, since we do not
   * need the blocks, us a standard MPI vector.
   */
  dealii::PETScWrappers::MPI::Vector             system_rhs;

  /**
   * The solution. This has a block structure but, since we do not
   * need the blocks, us a standard MPI vector. @note This vector has
   * hpghost elements.
   */
  dealii::PETScWrappers::MPI::Vector             solution;

  /**
   * An object to hold various run-time parameters that are specified
   * in an "prm file".
   */
  dealii::ParameterHandler parameters;

  /**
   * A table to handle a set of commonly interesting results.
   */
  dealii::TableHandler output_table;
};


/**
 * Constructor of the piezoelectric problem
 */
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
PiezoelectricProblem<dim, GroupSymm, ValueType>::PiezoelectricProblem ()
  :
  mpi_communicator (MPI_COMM_WORLD),

  pcout (std::cout, (dealii::Utilities::MPI::this_mpi_process (mpi_communicator) == 0)),
  
  triangulation (mpi_communicator,
		 typename dealii::Triangulation<dim>::MeshSmoothing
		 (dealii::Triangulation<dim>::smoothing_on_refinement |
		  dealii::Triangulation<dim>::smoothing_on_coarsening)),
  
  fe_q (dealii::FE_Q<dim> (2), dim, /* displacement       */
	dealii::FE_Q<dim> (1), 1),  /* electric potential */

  dof_handler (triangulation)
{}


/**
 * Destructor. This just ensures freeing some memory allocation.
 */
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
PiezoelectricProblem<dim, GroupSymm, ValueType>::~PiezoelectricProblem ()
{
  dof_handler.clear ();
}


/**
 * Obtain a complete list of the coefficients needed for this
 * calculation. First comes a declaration of the entries expected to
 * be find in the parameter file and then they are read into the
 * object parameters.
 */
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
PiezoelectricProblem<dim, GroupSymm, ValueType>::get_parameters ()
{
  // First declare the parameters that are expected to be found.
  parameters.declare_entry ("First-order elastic coefficients",
			    "-0.230",
			    dealii::Patterns::List (dealii::Patterns::Double (), 1),
			    "A list of the first-order elastic coefficients. " 
			    "Default is zinc-blende GaAs. ");

  parameters.declare_entry ("First-order dielectric coefficients",
			    "-0.230",
			    dealii::Patterns::List (dealii::Patterns::Double (), 1),
			    "A list of the first-order piezoelectric coefficients. " 
			    "Default is zinc-blende GaAs. ");

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

  parameters.declare_entry ("First-order spontaneous polarization coefficients",
			    "0.",
			    dealii::Patterns::List (dealii::Patterns::Double (), 1),
			    "A list of the first-order spontaneous polarization coefficients. " 
			    "Default is zinc-blende GaAs. ");

  parameters.declare_entry ("Bravais lattice dimensions",
   			    "1., 1., 1.",
   			    dealii::Patterns::List (dealii::Patterns::Double (), 1),
   			    "A list of the Bravais lattice dimensions. ");

  parameters.read_input (command_line.get_prm_file ());

  first_order_elastic_coefficients = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("First-order elastic coefficients")));

  first_order_dielectric_coefficients = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("First-order dielectric coefficients")));

  first_order_piezoelectric_coefficients = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("First-order piezoelectric coefficients")));

  second_order_piezoelectric_coefficients = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("Second-order piezoelectric coefficients")));

  first_order_spontaneous_polarization_coefficients = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("First-order spontaneous polarization coefficients")));

  bravais_lattice = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("Bravais lattice dimensions")));

}


/**
 * Generate the coarse grid that will be used for this calculation.
 */
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
PiezoelectricProblem<dim, GroupSymm, ValueType>::make_coarse_grid (const unsigned int n_refinement_cycles)
{
  dealii::GridGenerator::hyper_rectangle (triangulation,
					  dealii::Point<dim> (-20,-20,-20),
					  dealii::Point<dim> ( 20, 20, 20),
					  true);

  triangulation.refine_global (n_refinement_cycles);
}


/**
 * This function distribute coefficients onto the tensors of
 * coefficients. Before distribution, the tensors are reinitialised to
 * a zero state.
 */
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
PiezoelectricProblem<dim, GroupSymm, ValueType>::setup_coefficients ()
{
  // Initialise the first- 
  first_order_elastic_tensor.reinit ();
  first_order_dielectric_tensor.reinit ();
  first_order_piezoelectric_tensor.reinit ();

  // and second-order coefficient tensors,
  second_order_piezoelectric_tensor.reinit ();

  // and then distribute the first-
  first_order_elastic_tensor.distribute_coefficients (first_order_elastic_coefficients);
  first_order_dielectric_tensor.distribute_coefficients (first_order_dielectric_coefficients);
  first_order_piezoelectric_tensor.distribute_coefficients (first_order_piezoelectric_coefficients);

  // and second-order coefficients.
  second_order_piezoelectric_tensor.distribute_coefficients (second_order_piezoelectric_coefficients);

  pcout << "Tensors of coefficients:"
	<< std::endl
	<< "   First-order Elastic:             " << first_order_elastic_tensor
	<< std::endl
	<< "   First-order dielectric:          " << first_order_dielectric_tensor
	<< std::endl
	<< "   First-orfer piezoelectric:       " << first_order_piezoelectric_tensor
	<< std::endl
	<< "   First-orfer spont. polarization: " << first_order_spontaneous_polarization_tensor
	<< std::endl;
}


/**
 * Setup the linear algebra problem by initialising the boundary
 * constraints, matrices and vectors required for finding the
 * solution.
 */
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
PiezoelectricProblem<dim, GroupSymm, ValueType>::setup_system ()
{
  dof_handler.distribute_dofs (fe_q);

  locally_owned_dofs = dof_handler.locally_owned_dofs ();
  dealii::DoFTools::extract_locally_relevant_dofs (dof_handler,
						   locally_relevant_dofs);

  constraints.clear ();
  constraints.reinit (locally_relevant_dofs);
  dealii::DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  dealii::DoFTools::make_zero_boundary_constraints (dof_handler, constraints);
  constraints.close ();

  dealii::CompressedSimpleSparsityPattern csp (locally_relevant_dofs);
  dealii::DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
  dealii::SparsityTools::distribute_sparsity_pattern (csp,
						      dof_handler.n_locally_owned_dofs_per_processor (),
						      mpi_communicator,
						      locally_relevant_dofs);

  // @todo The system_matrix initialisation screws up with petsc-3.5.x
  system_matrix.reinit (locally_owned_dofs,
			locally_owned_dofs,
			csp,
			mpi_communicator);
  
  system_rhs.reinit (locally_owned_dofs, mpi_communicator);
  
  // The solution vector has ghost elements
  solution.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

}


/**
 * This function assembles the system matrix and rhight-hand-side
 * vector.
 */
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
PiezoelectricProblem<dim, GroupSymm, ValueType>::assemble_system ()
{
 
  dealii::QGauss<dim> face_quadrature_formula ();
  
}



/**
 * This is the run function, which wraps all of the above into a
 * single logical routine.
 */
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
PiezoelectricProblem<dim, GroupSymm, ValueType>::run ()
{

  // First find the parameters need for this calculation
  get_parameters ();
 
  // Then create the coarse grid
  make_coarse_grid (3);

  // Make coefficient tensors.
  setup_coefficients ();

  // Here comes the adaptive cycles
  for (unsigned int cycle=0; cycle<1; ++cycle)
    {
      
      pcout << "Grid:"
	    << std::endl
	    << "   Number of active cells:          "
	    << triangulation.n_global_active_cells ()
	    << " (on "
	    << triangulation.n_levels ()
	    << " levels)"
	    << std::endl;
      
      setup_system ();
      
      pcout << "   Number of degrees of freedom:    "
	    << dof_handler.n_dofs ()
	    << std::endl;


    }

}


int main (int argc, char **argv)
{

  // Initialise MPI as we allways do.
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, dealii::numbers::invalid_unsigned_int);

  try
    {
      dealii::deallog.depth_console (0);

      // Initialise the problem, 3d wurtzite, parse the command
      // line options and then run the problem.
      PiezoelectricProblem<3, nil::GroupSymmetry::Wurtzite> piezoelectric_problem;
      piezoelectric_problem.command_line.parse_command_line (argc, argv); 
      piezoelectric_problem.run ();
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
