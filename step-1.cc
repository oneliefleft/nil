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
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>

// Library-based headers.
#include "include/nil/group_symmetry.h"
#include "include/nil/dielectric_tensor.h"
#include "include/nil/elastic_tensor.h"
#include "include/nil/piezoelectric_tensor.h"
#include "include/nil/spontaneous_polarization_tensor.h"

#include "geometry_description.h"
#include "parameter_reader.h"

#include "piezoelectric_coefficients.h"

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
  PiezoelectricProblem (const std::string &parameter_file = "gaas.prm");

  /**
   * Destructor.
   */
  ~PiezoelectricProblem ();
  
  /**
   * Run the problem (a basic useage example)
   */
  void run ();
  
private:
  
  // First is a list of functions that belong to this class (they are
  // explained later on).
  void get_parameters ();

  void make_coarse_grid (const unsigned int n_refinement_cycles = 0);
  void make_boundary_constraints ();

  void setup_system ();
  void make_coefficient_tensors ();
  void assemble_system ();

  unsigned int solve ();

  void output_results (const unsigned int cycle) const;

  /**
   * A local copy of the MPI communicator.
   */
  MPI_Comm mpi_communicator;

  // Following that we have a list of the tensors that will be used in
  // this calculation. They are, first- 
  std::vector<nil::ElasticTensor<GroupSymm, 1, ValueType> >                 first_order_elastic_tensor;
  std::vector<nil::DielectricTensor<GroupSymm, 1, ValueType> >              first_order_dielectric_tensor;
  std::vector<nil::PiezoelectricTensor<GroupSymm, 1, ValueType> >           first_order_piezoelectric_tensor;
  std::vector<nil::SpontaneousPolarizationTensor<GroupSymm, 1, ValueType> > first_order_spontaneous_polarization_tensor;

  // and second-order piezoelectric tensors
  std::vector<nil::PiezoelectricTensor<GroupSymm, 2, ValueType> > second_order_piezoelectric_tensor;
 
  // Additionally, lists of coefficients are needed for first-order
  std::vector<std::vector<ValueType> > first_order_elastic_coefficients;
  std::vector<std::vector<ValueType> > first_order_dielectric_coefficients;
  std::vector<std::vector<ValueType> > first_order_piezoelectric_coefficients;
  std::vector<std::vector<ValueType> > first_order_spontaneous_polarization_coefficients;

  // and second-order tensors.
  std::vector<std::vector<ValueType> > second_order_piezoelectric_coefficients;

  /**
   * Sets of piezoelectric coefficients that are to be used in this
   * application - one for AlN and one for GaN.
   */
  nil::PiezoelectricCoefficients<nil::GroupSymmetry::Wurtzite, ValueType> aln_coefficients;
  nil::PiezoelectricCoefficients<nil::GroupSymmetry::Wurtzite, ValueType> gan_coefficients;

  /**
   * Mismatch strain tensor.
   */
  std::vector<dealii::Tensor<2, dim, ValueType> > lattice_mismatch_tensor;
  std::vector<std::vector<ValueType> >            lattice_coefficients;

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
  dealii::ParameterHandler prm_handler;

  /**
   * A table to handle a set of commonly interesting results.
   */
  dealii::TableHandler output_table;

  /**
   * The filename of the parameter file where material parameters are
   * kept. @note If no parameter file is given, the default GaAs is
   * used.
   */
  const std::string prm_file;

  /**
   * A separate (sub-)class that handles Postprocessing.
   */
  class Postprocessor;

  // A dummy number that counts how many ids we have.
  const unsigned int n_material_ids = 2;
};


/**
 * Constructor of the piezoelectric problem
 */
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
PiezoelectricProblem<dim, GroupSymm, ValueType>::PiezoelectricProblem (const std::string &parameter_file)
  :
  mpi_communicator (MPI_COMM_WORLD),

  pcout (std::cout, (dealii::Utilities::MPI::this_mpi_process (mpi_communicator) == 0)),
  
  triangulation (mpi_communicator,
		 typename dealii::Triangulation<dim>::MeshSmoothing
		 (dealii::Triangulation<dim>::smoothing_on_refinement |
		  dealii::Triangulation<dim>::smoothing_on_coarsening)),
  
  fe_q (dealii::FE_Q<dim> (2), dim, /* displacement       */
	dealii::FE_Q<dim> (1), 1),  /* electric potential */

  dof_handler (triangulation),

  prm_file (parameter_file)
{}


/**
 * Destructor. This just frees some memory allocation.
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

  nil::ParameterReader prm_reader (prm_handler, prm_file);
  prm_reader.read_parameters ();

#ifdef MANUAL_PARAMETR_HANDLER
  // First declare the parameters that are expected to be found.
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

  prm_handler.read_input (prm_file);
#endif

  for (unsigned int i=0; i<n_material_ids; ++i)
    {
      std::string subsection 
	= "Material id " + dealii::Utilities::int_to_string (i);
      
      prm_handler.enter_subsection (subsection);
      {
    
	first_order_elastic_coefficients.push_back (
	  dealii::Utilities::string_to_double
	  (dealii::Utilities::split_string_list (prm_handler.get ("First-order elastic coefficients"))));
	
	first_order_dielectric_coefficients.push_back (
	  dealii::Utilities::string_to_double
	  (dealii::Utilities::split_string_list (prm_handler.get ("First-order dielectric coefficients"))));
	
	first_order_piezoelectric_coefficients.push_back (
	  dealii::Utilities::string_to_double
	  (dealii::Utilities::split_string_list (prm_handler.get ("First-order piezoelectric coefficients"))));
	
	second_order_piezoelectric_coefficients.push_back ( 
	  dealii::Utilities::string_to_double
	  (dealii::Utilities::split_string_list (prm_handler.get ("Second-order piezoelectric coefficients"))));
	
	first_order_spontaneous_polarization_coefficients.push_back (
	  dealii::Utilities::string_to_double
	  (dealii::Utilities::split_string_list (prm_handler.get ("First-order spontaneous polarization coefficients"))));
	
	lattice_coefficients.push_back ( 
	  dealii::Utilities::string_to_double
	  (dealii::Utilities::split_string_list (prm_handler.get ("Bravais lattice coefficients"))));
      }
      prm_handler.leave_subsection ();
    }
}

  
/**
 * Postprocess piezoelectric solution.
 */
template <int dim, enum nil::GroupSymmetry group_symmetry = nil::None, typename ValueType = double>
class PiezoelectricProblem<dim, group_symmetry, ValueType>::Postprocessor 
  : 
  public dealii::DataPostprocessor<dim>
{
public:

  /**
   * Constructor.
   */  
  Postprocessor (const unsigned int partition);


  /**
   * Main postprocessing function. 
   */
  virtual
  void
  compute_derived_quantities_vector (const std::vector<dealii::Vector<ValueType> >                       &uh,
				     const std::vector<std::vector<dealii::Tensor<1, dim, ValueType> > > &duh,
				     const std::vector<std::vector<dealii::Tensor<2, dim, ValueType> > > &dduh,
				     const std::vector<dealii::Point<dim, ValueType> >                   &normals,
				     const std::vector<dealii::Point<dim, ValueType> >                   &evaluation_points,
				     std::vector<dealii::Vector<ValueType> >                             &computed_quantities) const;

  
  /**
   * Get an the names of each of the components of the output vector.
   */  
  virtual 
  std::vector<std::string> get_names () const;


  /**
   * Get an interpretation of what the components of the output vector mean.
   */  
  virtual
  std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation> 
  get_data_component_interpretation () const;
  

  /**
   * Get the update flags for the subdomain \f$h\f$ (ie. finite
   * element) thata re needed to postprocess the solution in
   * quadrature.
   */
  virtual 
  dealii::UpdateFlags get_needed_update_flags () const;
  
private:
  
  const unsigned int partition;

};


template <int dim, enum nil::GroupSymmetry group_symmetry, typename ValueType>
PiezoelectricProblem<dim, group_symmetry, ValueType>::Postprocessor::Postprocessor (const unsigned int partition)
  :
  partition (partition)
{}


template <int dim, enum nil::GroupSymmetry group_symmetry, typename ValueType>
std::vector<std::string>
PiezoelectricProblem<dim, group_symmetry, ValueType>::Postprocessor::get_names () const
{
  std::vector<std::string> solution_names;

  // First deal with the displacement vector (u_i).  
  for (unsigned int d=0; d<dim; ++d)
    {
      const std::string u_name = "u_" + dealii::Utilities::int_to_string (d);
      solution_names.push_back (u_name);
    }

  // and then the potential (phi).
  solution_names.push_back ("phi");

  // solution_names.push_back ("partition");
  
  return solution_names;
}


template <int dim, enum nil::GroupSymmetry group_symmetry, typename ValueType>
std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
PiezoelectricProblem<dim, group_symmetry, ValueType>::Postprocessor::get_data_component_interpretation () const
{
  // First deal with the displacement vector (u_i).
  std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
    interpretation (dim, dealii::DataComponentInterpretation::component_is_part_of_vector);
  
  // and then the potential (phi).
  interpretation.push_back (dealii::DataComponentInterpretation::component_is_scalar);
  
  return interpretation;
}


template <int dim, enum nil::GroupSymmetry group_symmetry, typename ValueType>
dealii::UpdateFlags
PiezoelectricProblem<dim, group_symmetry, ValueType>::Postprocessor::get_needed_update_flags() const
{
  return 
    dealii::update_values    | 
    dealii::update_gradients | 
    dealii::update_q_points;
}


template <int dim, enum nil::GroupSymmetry group_symmetry, typename ValueType>
void
PiezoelectricProblem<dim, group_symmetry, ValueType>::Postprocessor::
compute_derived_quantities_vector (const std::vector<dealii::Vector<ValueType> >                       &uh,
				   const std::vector<std::vector<dealii::Tensor<1, dim, ValueType> > > &duh,
				   const std::vector<std::vector<dealii::Tensor<2, dim, ValueType> > > &/* dduh              */,
				   const std::vector<dealii::Point<dim, ValueType> >                   &/* normals           */,
				   const std::vector<dealii::Point<dim, ValueType> >                   &/* evaluation_points */,
				   std::vector<dealii::Vector<ValueType> >                             &computed_quantities) const
{
  const unsigned int n_quadrature_points = uh.size ();

  Assert (duh.size() == n_quadrature_points, dealii::ExcInternalError ());
  Assert (computed_quantities.size () == n_quadrature_points, dealii::ExcInternalError ());
  Assert (uh[0].size() == dim+1, dealii::ExcInternalError());

  for (unsigned int q=0; q<n_quadrature_points; ++q)
     {
       // First deal with the displacement vector (u_i).  
       for (unsigned int d=0; d<dim; ++d)
	 computed_quantities[q](d) = uh[q](d);

       // and then the potential (phi).
       computed_quantities[q](dim) = uh[q](dim+1);

       //     Tensor<2,dim> grad_u;
       //     for (unsigned int d=0; d<dim; ++d)
       // 	grad_u[d] = duh[q][d];
       //     const SymmetricTensor<2,dim> strain_rate = symmetrize (grad_u);
       //     computed_quantities[q](dim+2) = 2 * EquationData::eta *
       // 	strain_rate * strain_rate;
       //     computed_quantities[q](dim+3) = partition;
     }
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
PiezoelectricProblem<dim, GroupSymm, ValueType>::make_coefficient_tensors ()
{
  // Initialise the first- 
  first_order_elastic_tensor.resize (n_material_ids);
  first_order_dielectric_tensor.resize (n_material_ids);
  first_order_piezoelectric_tensor.resize (n_material_ids);
  first_order_spontaneous_polarization_tensor.resize (n_material_ids);

  // and second-order coefficient tensors,
  second_order_piezoelectric_tensor.resize (n_material_ids);

  for (unsigned int i=0; i<n_material_ids; ++i)
    {
      // and then distribute the first-
      first_order_elastic_tensor[i].distribute_coefficients (first_order_elastic_coefficients[i]);
      first_order_dielectric_tensor[i].distribute_coefficients (first_order_dielectric_coefficients[i]);
      first_order_piezoelectric_tensor[i].distribute_coefficients (first_order_piezoelectric_coefficients[i]);
      first_order_spontaneous_polarization_tensor[i].distribute_coefficients (first_order_spontaneous_polarization_coefficients[i]);
      
      // and second-order coefficients.
      second_order_piezoelectric_tensor[i].distribute_coefficients (second_order_piezoelectric_coefficients[i]);
    }

  // mismatch strain
  lattice_mismatch_tensor.resize (n_material_ids);
  lattice_mismatch_tensor[0].clear ();
  lattice_mismatch_tensor[1].clear ();

  const double delta_a = lattice_coefficients[1][0]/lattice_coefficients[0][0]-1.;
  const double delta_c = lattice_coefficients[1][1]/lattice_coefficients[0][1]-1.;

  lattice_mismatch_tensor[1][0][0] = lattice_mismatch_tensor[1][1][1] = delta_a;
  lattice_mismatch_tensor[1][2][2] = delta_c;
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
 * This function assembles the system matrix and right-hand-side
 * vector.
 */
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
void 
PiezoelectricProblem<dim, GroupSymm, ValueType>::assemble_system ()
{
  // @todo: The number of quadrature points should be decided
  // elsewhere. @note: overspecify the number of quadrature points.
  dealii::QGauss<dim> quadrature_formula (9);
  
  dealii::FEValues<dim> fe_values (fe_q, quadrature_formula,
				   dealii::update_quadrature_points | 
				   dealii::update_values            | 
				   dealii::update_gradients         | 
				   dealii::update_JxW_values);
  
  const unsigned int n_dofs_per_cell = fe_q.dofs_per_cell;
  const unsigned int n_q_points      = quadrature_formula.size ();
  
  dealii::FullMatrix<ValueType> cell_matrix (n_dofs_per_cell, n_dofs_per_cell);
  dealii::Vector<ValueType>     cell_rhs    (n_dofs_per_cell);

  std::vector<dealii::types::global_dof_index> local_dof_indices (n_dofs_per_cell);

  // Get a copy of the geometry for the quantum dot
  nil::GeometryDescription::HyperCube<dim, ValueType> hyper_cube (-5,5);
  std::vector<ValueType> cell_hyper_cube (n_q_points);

  // So-called "finite elemnt views" that define parts of the finite
  // elemnt system connected with displacements (u on position
  // [0,1,2]) and electric potential (phi on position 3).
  const dealii::FEValuesExtractors::Vector u   (0);
  const dealii::FEValuesExtractors::Scalar phi (dim);

  typename dealii::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  
  for (; cell != endc; ++cell)
    if (cell->is_locally_owned ())
      {

	// Initialise objects that are needed to be in a clean state
	// on this cell.
	fe_values.reinit (cell);
	cell_matrix = 0;
	cell_rhs    = 0;

	// Obtain the material pattern on quandrature points.
	hyper_cube.value_list (fe_values.get_quadrature_points (), cell_hyper_cube);

	for (unsigned int q_point = 0; q_point<n_q_points; ++q_point)
	  {

	    // This is the material id 1 in the dot and 0 otherwise.
	    const unsigned int material_id = cell_hyper_cube[q_point];

	    for (unsigned int i=0; i<n_dofs_per_cell; ++i)
	      {

		// obtain symmetric gradient for the ith q_point.
                const dealii::Tensor<2, dim> u_i_grad   = fe_values[u].symmetric_gradient (i, q_point);
                const dealii::Tensor<1, dim> phi_i_grad = fe_values[phi].gradient (i, q_point);

		for (unsigned int j=0; j<n_dofs_per_cell; ++j)
		  {
		    
		    // obtain symmetric gradient for the jth q_point.
		    const dealii::Tensor<2, dim> u_j_grad   = fe_values[u].symmetric_gradient (j, q_point);
		    const dealii::Tensor<1, dim> phi_j_grad = fe_values[phi].gradient (j, q_point);
		    
		    cell_matrix (i,j) +=
                       (contract (u_i_grad, first_order_elastic_tensor[material_id], 
				  u_j_grad)
			+
			contract (phi_i_grad, first_order_dielectric_tensor[material_id], 
				  phi_j_grad))
		      *
		      fe_values.JxW (q_point);
		    
		    
		  } // dof j
		
		// @todo By commenting out the first line, the solver
		// converges (no -nan).
		cell_rhs (i) -= 
		  (contract (u_i_grad, first_order_elastic_tensor[material_id], 
		   	     lattice_mismatch_tensor[material_id])
		   +
		   contract (phi_i_grad, first_order_spontaneous_polarization_tensor[material_id]))
		  
		  *
		  fe_values.JxW (q_point);
		
	      } // dof i

	  } // q_point

	cell->get_dof_indices(local_dof_indices);

	constraints.distribute_local_to_global (cell_matrix, local_dof_indices,
						system_matrix);

	constraints.distribute_local_to_global (cell_rhs, local_dof_indices,
						system_rhs);

      } // cell

  system_matrix.compress (dealii::VectorOperation::add);
  system_rhs.compress (dealii::VectorOperation::add);

}


/**
 * Solve the linear algebra system. @note In the general case, the
 * matrix is not symmetric.
 */
template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
unsigned int
PiezoelectricProblem<dim, GroupSymm, ValueType>::solve ()
{

  // Ultimately, a vector is required that stores not only the
  // elements of the solution for degrees of freedom the current
  // processor owns, but also all other locally relevant degrees of
  // freedom. The solver itself needs a vector that is uniquely split
  // between processors, without any overlap. Therefore create a
  // vector at the beginning of this function that has these
  // properties, use it to solve the linear system, and only assign it
  // to the vector we want at the very end. This last step ensures
  // that all ghost elements are also copied as necessary.
  dealii::PETScWrappers::MPI::Vector distributed_solution (locally_owned_dofs, mpi_communicator);
  dealii::SolverControl solver_control (system_matrix.m (), 1e-09*system_rhs.l2_norm ());
  dealii::PETScWrappers::PreconditionBlockJacobi preconditioner (system_matrix);
  dealii::PETScWrappers::SolverBicgstab solver (solver_control, mpi_communicator);
  // dealii::PETScWrappers::SolverGMRES solver (solver_control, mpi_communicator);

  solver.solve (system_matrix, distributed_solution, system_rhs, 
		preconditioner);

  constraints.distribute (distributed_solution);
  solution = distributed_solution;

  return solver_control.last_step ();
}


template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
void
PiezoelectricProblem<dim, GroupSymm, ValueType>::output_results (const unsigned int cycle) const
{
  Postprocessor postprocessor (dealii::Utilities::MPI::this_mpi_process (mpi_communicator));

  dealii::DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, postprocessor);
  data_out.build_patches ();

  const std::string filename 
    = ("solution-" +
       dealii::Utilities::int_to_string (cycle, 4) +
       "." +
       dealii::Utilities::int_to_string (triangulation.locally_owned_subdomain(), 4) +
       ".vtu");
  
  std::ofstream output (filename.c_str());

  data_out.write_vtu (output);

  // For calculations in parallel it is convenient to have a master
  // record written. This allows paravew (and friends) to read the
  // entire structure by binding blocks (*.vtu) into a single
  // structure (*.visit). This only needs to be done by one processor.
  if (dealii::Utilities::MPI::this_mpi_process (mpi_communicator) == 0)
    {
      std::vector<std::string> filenames;

      // The filenames to bind into the master record should match the
      // filenames used above.
      for (unsigned int i=0; i<dealii::Utilities::MPI::n_mpi_processes (mpi_communicator); ++i)
	{
	  filenames.push_back ("solution-" +
			       dealii::Utilities::int_to_string (cycle, 4) +
			       "." +
			       dealii::Utilities::int_to_string (i, 4) +
			       ".vtu");
	}
      pcout << std::endl;

      const std::string
	pvtu_master_filename = ("solution-" +
				dealii::Utilities::int_to_string (cycle, 4) +
				".pvtu");
  
      std::ofstream pvtu_master (pvtu_master_filename.c_str ());
      data_out.write_pvtu_record (pvtu_master, filenames);

      // Create master record for visit
      const std::string visit_master_filename 
	= ("solution-" +
	   dealii::Utilities::int_to_string (cycle, 4) +
	   ".visit");

      std::ofstream visit_master (visit_master_filename.c_str ());
      data_out.write_visit_record (visit_master, filenames);
    }

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

  // Make coefficient tensors.
  make_coefficient_tensors ();

  for (unsigned int i=0; i<n_material_ids; ++i)
    { 
      pcout << std::endl
	    << "Tensors of coefficients:"
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
	    << "   1st-order spont. polarization:   " << first_order_spontaneous_polarization_tensor[i]
	    << std::endl
	    << "   2nd-order piezoelectric:         " << second_order_piezoelectric_tensor[i]
	    << std::endl;
    }

  // Then create the coarse grid
  make_coarse_grid (2);

  // Here comes the adaptive cycles
  for (unsigned int cycle=0; cycle<1; ++cycle)
    {
      
      pcout << std::endl
	    << "Grid:"
	    << std::endl
	    << "   Number of active cells:          "
	    << triangulation.n_global_active_cells ()
	    << " (on "
	    << triangulation.n_levels ()
	    << " levels)"
	    << std::endl;
      
      setup_system ();
      
      pcout << std::endl
	    << "System:"
	    << std::endl
	    << "   Number of degrees of freedom:    "
	    << dof_handler.n_dofs ()
	    << std::endl;

      assemble_system ();

      pcout << std::endl
	    << "Linear algebra system:"
	    << std::endl
	    << "   Number of non-zero elements:     "
	    << system_matrix.n_nonzero_elements ()
	    << std::endl
	    << "   System rhs l2-norm:              "
	    << system_rhs.l2_norm ()
	    << std::endl << std::flush;

      // std::cout << std::endl
      // 		<< "System matrix:"
      // 		<< std::endl;
      // system_matrix.print (std::cout);
      // std::cout << std::endl
      // 		<< "System rhs:"
      // 		<< std::endl;
      // system_rhs.print (std::cout);
      // std::cout << std::endl;

      const unsigned int n_iterations = solve ();

      pcout << "Solver converged in:                "
	    << n_iterations << " iterations"
	    << std::endl;

      output_results (cycle);
    }

}


int main (int argc, char **argv)
{

  // Initialise MPI as we allways do.
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, dealii::numbers::invalid_unsigned_int);

  try
    {
      dealii::deallog.depth_console (0);

      // Initialise the problem, 3d wurtzite, (default double).
      PiezoelectricProblem<3, nil::GroupSymmetry::Wurtzite> piezoelectric_problem ("step-1.prm");
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
