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

#ifndef __nil_piezoelectric_model_h
#define __nil_piezoelectric_model_h

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

namespace nil
{

namespace Piezoelectric
{

  /**
   * This class can setup, solve, and output the results of the
   * piezoelectric problem.
   *
   * @author Toby D. Young 2012, 2014
   */ 
  template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType = double>
    class Model
    {
    public:
    
    /**
     * Constructor. Take a parameter file name if any.
     */
    Model (const std::string &parameter_file = "");
    
    /**
     * Destructor.
     */
    ~Model ();
    
    /**
     * Run the problem.
     */
    void run ();
    
  private:
    
    // First is a list of functions that belong to this class (they are
    // explained later on).
    void get_parameters ();
    
    void make_coarse_grid (const unsigned int n_refinement_cycles = 0);
    void make_boundary_constraints ();
    
    void setup_system ();
    void setup_coefficient_tensors ();
    void assemble_system ();
    
    unsigned int solve ();
    void refine_grid ();
    
    void output_results (const unsigned int cycle) const;
    void output_material_id (const unsigned int cycle) const;
    
    // A local copy of the MPI communicator.
    MPI_Comm mpi_communicator;
    
    // Following that we have a list of the tensors that will be used
    // in this calculation. They are, first-
    std::vector<nil::ElasticTensor<GroupSymm, 1, ValueType> >       first_order_elastic_tensor;
    std::vector<nil::DielectricTensor<GroupSymm, 1, ValueType> >    first_order_dielectric_tensor;
    std::vector<nil::PiezoelectricTensor<GroupSymm, 1, ValueType> > first_order_piezoelectric_tensor;
    std::vector<nil::PolarelectricTensor<GroupSymm, 1, ValueType> > first_order_polarelectric_tensor;
    
    // and second-order tensors of coefficients
    std::vector<nil::PiezoelectricTensor<GroupSymm, 2, ValueType> > second_order_piezoelectric_tensor;
    
    // Additionally, lists of coefficients are needed for first-order
    std::vector<std::vector<ValueType> > first_order_elastic_coefficients;
    std::vector<std::vector<ValueType> > first_order_dielectric_coefficients;
    std::vector<std::vector<ValueType> > first_order_piezoelectric_coefficients;
    std::vector<std::vector<ValueType> > first_order_polarelectric_coefficients;
    
    // and second-order tensors.
    std::vector<std::vector<ValueType> > second_order_piezoelectric_coefficients;
    
    // Mismatch strain tensor.
    std::vector<dealii::Tensor<2, dim, ValueType> > lattice_mismatch_tensor;
    std::vector<std::vector<ValueType> >            lattice_coefficients;
    
    
    // A I<code>deal.II</code> hack that outputs to the first processor
    // only (useful for output in parallel calculations).
    dealii::ConditionalOStream pcout;
    
    // A parallel distributed triangulation.
    dealii::parallel::distributed::Triangulation<dim> triangulation;
    
    // The finite element and linear algebra system.
    const dealii::FESystem<dim>              fe_q;
    dealii::DoFHandler<dim>                  dof_handler;
    dealii::ConstraintMatrix                 constraints;
    dealii::IndexSet                         locally_owned_dofs;
    dealii::IndexSet                         locally_relevant_dofs;
    
    // Objects for linear algebra calculation
#ifdef USE_PETSC
    dealii::PETScWrappers::MPI::SparseMatrix system_matrix;
    dealii::PETScWrappers::MPI::Vector       system_rhs;
    dealii::PETScWrappers::MPI::Vector       solution;
#else
    dealii::TrilinosWrappers::SparseMatrix   system_matrix;
    dealii::TrilinosWrappers::MPI::Vector    system_rhs;
    dealii::TrilinosWrappers::MPI::Vector    solution;
#endif
    

    // An object to hold various run-time parameters that are specified
    // in a "prm file".
    dealii::ParameterHandler prm_handler;
    const std::string        prm_file;
    
    // Piezoelectric postprocessor.
    class Postprocessor;
    
    // A dummy number that counts how many material ids we have.
    const unsigned int n_material_ids;
  };
  
} // namespace Piezoelectric

} // namespace nil

#endif // __nil_piezoelectric_model_h
