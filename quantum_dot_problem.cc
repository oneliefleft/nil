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
#include <deal.II/lac/sparsity_tools.h>

// USE_PETSC
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_vector.h>

// USE_TRILINOS
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

// USE_SLEPC
#include <deal.II/lac/slepc_solver.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

// Library-based headers.
#include "include/nil/base/group_symmetry.h"

#include "include/nil/dielectric_tensor.h"
#include "include/nil/elastic_tensor.h"
#include "include/nil/piezoelectric_tensor.h"
#include "include/nil/polarelectric_tensor.h"

#include "include/nil/geometry_function.h"

#include "include/nil/parameter_reader.h"

#include "piezoelectric_model.h"
#include "piezoelectric_coefficients.h"
#include "valance_band_model.h"

#include <fstream>
#include <iostream>


template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
class QuantumDotProblem
{
public:

  /**
   * Constructor.
   */  
  QuantumDotProblem ();

  /**
   * Destructor.
   */
  ~QuantumDotProblem ();

  void run ();

  /**
   * The number of material ids.
   */
  unsigned int n_material_ids;

private:

  /**
   * A local copy of the MPI communicator.
   */
  MPI_Comm mpi_communicator;

  /**
   * A I<code>deal.II</code> hack that outputs to the first processor
   * only (useful for output in parallel calculations).
   */
  dealii::ConditionalOStream pcout;  

  // A parallel distributed triangulation.
  dealii::parallel::distributed::Triangulation<dim> triangulation;
};


template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
QuantumDotProblem<dim, GroupSymm, ValueType>::QuantumDotProblem ()
  :
  mpi_communicator (MPI_COMM_WORLD),

  pcout (std::cout, (dealii::Utilities::MPI::this_mpi_process (mpi_communicator) == 0)),
  
  triangulation (mpi_communicator,
		 typename dealii::Triangulation<dim>::MeshSmoothing
		 (dealii::Triangulation<dim>::smoothing_on_refinement |
		  dealii::Triangulation<dim>::smoothing_on_coarsening))
{
  n_material_ids = 2;
}


template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
QuantumDotProblem<dim, GroupSymm, ValueType>::~QuantumDotProblem ()
{}


template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
void QuantumDotProblem<dim, GroupSymm, ValueType>::run ()
{
  // Call the constructor for the piezoelectric problem.
  nil::Piezoelectric::Model<dim, GroupSymm, ValueType> piezoelectric_model 
    (triangulation,
     mpi_communicator);

  // First find the parameters need for this calculation
  piezoelectric_model.get_parameters ("piezoelectric.prm");
  
  // Make coefficient tensors.
  piezoelectric_model.setup_coefficient_tensors ();
  
#ifdef VERBOSE
  for (unsigned int i=0; i<n_material_ids; ++i)
    { 
      pcout << std::endl
	    << "Tensors of coefficients for material " << i << ":"
	    << std::endl
	    << "   Lattice coefficients:            ";
      for (unsigned int j=0; j<lattice_coefficients[i].size (); ++j)
	pcout << lattice_coefficients[i][j] << " ";
      pcout << std::endl
	    << "   Lattice mismatch:                " << lattice_mismatch_tensor[i]
	    << std::endl
	    << "   1st-order elastic:               " << first_order_elastic_tensor[i]
	    << std::endl
	    << "   1st-order dielectric:            " << first_order_dielectric_tensor[i]
	    << std::endl
	    << "   1st-order piezoelectric:         " << first_order_piezoelectric_tensor[i]
	    << std::endl
	    << "   1st-order spont. polarization:   " << first_order_polarelectric_tensor[i]
	    << std::endl
	    << "   2nd-order piezoelectric:         " << second_order_piezoelectric_tensor[i]
	    << std::endl;
    }
#endif
  
  // Then create the coarse grid
  piezoelectric_model.make_coarse_grid (4);
  
  // Here comes the adaptive cycles
  const unsigned int n_cycles = 4;
      
  for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
    {
      
      pcout << std::endl
	    << "Cycle:                              "
	    << cycle
	    << std::endl;

#ifdef VERBOSE      
      pcout << "   Number of active cells:          "
	    << triangulation.n_global_active_cells ()
	    << " (on "
	    << triangulation.n_levels ()
	    << " levels)"
	    << std::endl;
#endif
      
      // Distribute degrees of freedom and reinitialize matrices and
      // vectors.
      piezoelectric_model.setup_system ();
      
#ifdef VERBOSE
      pcout << "   Number of degrees of freedom:    "
	    << dof_handler.n_dofs ()
	    << std::endl;
#endif      

      // Assemble the matrix and rhs vector.
      piezoelectric_model.assemble_system ();
      
#ifdef VERBOSE
      pcout << "   Number of non-zero elements:     "
	    << system_matrix.n_nonzero_elements ()
	    << std::endl
	    << "   System rhs l2-norm:              "
	    << system_rhs.l2_norm ()
	    << std::endl;
#endif 
      
      // Solve the problem.
      const unsigned int n_iterations = piezoelectric_model.solve ();
      
      pcout << "   Solver converged in:             "
	    << n_iterations << " iterations"
	    << std::endl;
      
      // Write derived quantities to files.
      piezoelectric_model.output_results (cycle);
      piezoelectric_model.output_material_id (cycle);
      
      // Finally use Kelly's error estimate to refine the grid.
      piezoelectric_model.refine_grid ();
    }

  // Initialise the model, 3d wurtzite, (default double).
  nil::ValanceBand::Model<3, nil::GroupSymmetry::Wurtzite> six_band_hole_model 
  (triangulation,
   mpi_communicator);

}


  
int main (int argc, char **argv)
{

  // Initialise MPI as we allways do.
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, dealii::numbers::invalid_unsigned_int);

  try
    {
      dealii::deallog.depth_console (0);

      // Initialise the model, 3d wurtzite, (default double).
      QuantumDotProblem<3, nil::GroupSymmetry::Wurtzite, double> quantum_dot_problem;
      quantum_dot_problem.run ();
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
