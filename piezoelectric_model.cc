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

#include "piezoelectric_model.h"

namespace nil
{
  
  namespace Piezoelectric
  {
    
    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    Model<dim, GroupSymm, ValueType>::Model (const std::string &parameter_file)
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
      
      prm_file (parameter_file),
      
      n_material_ids (2)
    {}
    
    
    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    Model<dim, GroupSymm, ValueType>::~Model ()
    {
      // free some meory
      dof_handler.clear ();
    }
    
    
    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    void 
    Model<dim, GroupSymm, ValueType>::get_parameters ()
    {
      // The parameter file should not be an empty field if this
      // function is called.
      AssertThrow (prm_file!="",
		   dealii::ExcMessage ("The parameter file name is empty or does not exist"));
      
      // Assign parameters to the parameter handler and check sanity
      // of t he parameter. @note If no file exists an unusal default
      // is made.
      nil::ParameterReader prm_reader (prm_handler, prm_file);
      prm_reader.read_parameters ();
      
      // then read in the parameters and store them
      for (unsigned int i=0; i<n_material_ids; ++i)
	{
	  std::string subsection 
	    = "Material id " + dealii::Utilities::int_to_string (i);
	  
	  prm_handler.enter_subsection (subsection);
	  {
	    
	    first_order_elastic_coefficients.push_back 
	      (dealii::Utilities::string_to_double
	       (dealii::Utilities::split_string_list (prm_handler.get ("First-order elastic coefficients"))));
	    
	    first_order_dielectric_coefficients.push_back 
	      (dealii::Utilities::string_to_double
	       (dealii::Utilities::split_string_list (prm_handler.get ("First-order dielectric coefficients"))));
	    
	    first_order_piezoelectric_coefficients.push_back 
	      (dealii::Utilities::string_to_double
	       (dealii::Utilities::split_string_list (prm_handler.get ("First-order piezoelectric coefficients"))));
	    
	    second_order_piezoelectric_coefficients.push_back 
	      (dealii::Utilities::string_to_double
	       (dealii::Utilities::split_string_list (prm_handler.get ("Second-order piezoelectric coefficients"))));
	    
	    first_order_polarelectric_coefficients.push_back 
	      (dealii::Utilities::string_to_double
	       (dealii::Utilities::split_string_list (prm_handler.get ("First-order spontaneous polarization coefficients"))));
	    
	    lattice_coefficients.push_back 
	      (dealii::Utilities::string_to_double
	       (dealii::Utilities::split_string_list (prm_handler.get ("Bravais lattice coefficients"))));
	  }
	  prm_handler.leave_subsection ();
	}
    }
    
    
    /**
     * Postprocess piezoelectric solution.
     */
    template <int dim, enum nil::GroupSymmetry group_symmetry = nil::None, typename ValueType = double>
    class Model<dim, group_symmetry, ValueType>::Postprocessor 
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
       * Get an the names of each of the components of the output
       * vector.
       */  
      virtual 
      std::vector<std::string> get_names () const;
      
      /**
       * Get an interpretation of what the components of the output
       * vector mean.
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
    Model<dim, group_symmetry, ValueType>::Postprocessor::Postprocessor (const unsigned int partition)
      :
      partition (partition)
    {}
    
    
    template <int dim, enum nil::GroupSymmetry group_symmetry, typename ValueType>
    std::vector<std::string>
    Model<dim, group_symmetry, ValueType>::Postprocessor::get_names () const
    {
      std::vector<std::string> solution_names;
      
      // First deal with the displacement vector (u_i).  
      for (unsigned int d=0; d<dim; ++d)
	{
	  const std::string u_name = "u_" + dealii::Utilities::int_to_string (d+1);
	  solution_names.push_back (u_name);
	}
      
      // and then the potential (phi).
      solution_names.push_back ("phi");
      
      // solution_names.push_back ("partition");
      
    return solution_names;
    }
    
    
    template <int dim, enum nil::GroupSymmetry group_symmetry, typename ValueType>
    std::vector<dealii::DataComponentInterpretation::DataComponentInterpretation>
    Model<dim, group_symmetry, ValueType>::Postprocessor::get_data_component_interpretation () const
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
    Model<dim, group_symmetry, ValueType>::Postprocessor::get_needed_update_flags() const
    {
      return 
	dealii::update_quadrature_points | 
	dealii::update_values            | 
	dealii::update_gradients;
    }
    
    
    template <int dim, enum nil::GroupSymmetry group_symmetry, typename ValueType>
    void
    Model<dim, group_symmetry, ValueType>::Postprocessor::
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
	  computed_quantities[q](dim) = uh[q](dim);
	  
	  dealii::Tensor<2, dim, ValueType> grad_u;
	  
	  for (unsigned int d=0; d<dim; ++d)
	    grad_u[d] = duh[q][d];
	  
	  const dealii::SymmetricTensor<2, dim, ValueType> strain = symmetrize (grad_u);
	}
    }
    
        
    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    void 
    Model<dim, GroupSymm, ValueType>::make_coarse_grid (const unsigned int n_refinement_cycles)
    {
      dealii::GridGenerator::hyper_rectangle (triangulation,
					      dealii::Point<dim> (-20,-20,-20),
					      dealii::Point<dim> ( 20, 20, 20),
					      true);
      
      triangulation.refine_global (n_refinement_cycles);
    }
    
    
    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    void 
    Model<dim, GroupSymm, ValueType>::setup_coefficient_tensors ()
    {
      // Initialise the first- 
      first_order_elastic_tensor.resize (n_material_ids);
      first_order_dielectric_tensor.resize (n_material_ids);
      first_order_piezoelectric_tensor.resize (n_material_ids);
      first_order_polarelectric_tensor.resize (n_material_ids);
      
      // and second-order coefficient tensors,
      second_order_piezoelectric_tensor.resize (n_material_ids);
      
      for (unsigned int i=0; i<n_material_ids; ++i)
	{
	  // and then distribute the first-
	  first_order_elastic_tensor[i].distribute_coefficients (first_order_elastic_coefficients[i]);
	  first_order_dielectric_tensor[i].distribute_coefficients (first_order_dielectric_coefficients[i]);
	  first_order_piezoelectric_tensor[i].distribute_coefficients (first_order_piezoelectric_coefficients[i]);
	  first_order_polarelectric_tensor[i].distribute_coefficients (first_order_polarelectric_coefficients[i]);
	  
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
    
    
    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    void 
    Model<dim, GroupSymm, ValueType>::setup_system ()
    {
      dof_handler.distribute_dofs (fe_q);
      
      locally_owned_dofs = dof_handler.locally_owned_dofs ();
      dealii::DoFTools::extract_locally_relevant_dofs (dof_handler,
						       locally_relevant_dofs);
      
      // Deal with hanging-node and boundary constraints.
      constraints.clear ();
      constraints.reinit (locally_relevant_dofs);
      dealii::DoFTools::make_hanging_node_constraints (dof_handler, constraints);
      dealii::DoFTools::make_zero_boundary_constraints (dof_handler, constraints);
      constraints.close ();
      
      // Obtain a compressed sparsity pattern for sparse-matrix memory
      // allocation.
      dealii::CompressedSimpleSparsityPattern csp (locally_relevant_dofs);
      dealii::DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
      dealii::SparsityTools::distribute_sparsity_pattern (csp,
							  dof_handler.n_locally_owned_dofs_per_processor (),
							  mpi_communicator,
							  locally_relevant_dofs);
      
      // @todo The system_matrix initialisation screws up with petsc-3.5.x.
      system_matrix.reinit (locally_owned_dofs,
			    locally_owned_dofs,
			    csp,
			    mpi_communicator);
      
      system_rhs.reinit (locally_owned_dofs, mpi_communicator);
      
      // The solution vector has ghost elements.
      solution.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
      
    }
    
    
    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    void 
    Model<dim, GroupSymm, ValueType>::assemble_system ()
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
      nil::GeometryFunction::HalfHyperBall<dim, ValueType> geometry (5., dealii::Point<dim, ValueType> ());
      std::vector<ValueType> cell_geometry (n_q_points);
      
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
	    
	    // Initialise objects that are needed to be in a clean
	    // state on this cell.
	    fe_values.reinit (cell);
	    cell_matrix = 0;
	    cell_rhs    = 0;
	    
	    // Obtain the material pattern on quandrature points.
	    geometry.value_list (fe_values.get_quadrature_points (), cell_geometry);
	    
	    for (unsigned int q_point = 0; q_point<n_q_points; ++q_point)
	      {
	      
		// This is the material id 1 in the dot and 0 otherwise.
		const unsigned int material_id = cell_geometry[q_point];
		
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
			
			// The cell matrix is a block matrix of the form:
			//
			// | A B |
			// | C D |
			//
			// Where the components of A are the elastic
			// part, B are the reverse piezoelectric part, C
			// are the piezoelectric part, and D are the
			// dielectric part. In that order, the cell
			// matrix is assembled below.
			cell_matrix (i,j) +=
			  (contract (u_i_grad, first_order_elastic_tensor[material_id], 
				     u_j_grad)
			   -
			   contract (u_i_grad, first_order_piezoelectric_tensor[material_id], 
				     phi_j_grad)
			   +
			   contract (phi_i_grad, first_order_dielectric_tensor[material_id], 
				     phi_j_grad)
			   + 
			   contract (phi_i_grad, first_order_piezoelectric_tensor[material_id], 
				     u_j_grad))
			  
			  *
			  fe_values.JxW (q_point);
			
			
		      } // dof j
		    
		    
		    // The cell vector is a block vector of the form:
		    // 
		    // |x|
		    // |y|
		    //
		    // where the components x are the lattice mismatch
		    // strain and y are the polarelectric field. In that
		    // order, the cell vector is assembled below.
		    cell_rhs (i) += 
		      (contract (u_i_grad, first_order_elastic_tensor[material_id], 
				 lattice_mismatch_tensor[material_id])
		       +
		       contract (phi_i_grad, first_order_polarelectric_tensor[material_id]))
		      
		      *
		      fe_values.JxW (q_point);
		    
		  } // dof i
		
	      } // q_point
	    
	    cell->get_dof_indices(local_dof_indices);
	    
	    // Distribute the cell matrix onto the distributed system
	    // matrix and the cell vector onto the distributed rhs vector.
	    constraints.distribute_local_to_global (cell_matrix, local_dof_indices,
						    system_matrix);
	    
	    constraints.distribute_local_to_global (cell_rhs, local_dof_indices,
						    system_rhs);
	    
	  } // cell
      
      system_matrix.compress (dealii::VectorOperation::add);
      system_rhs.compress (dealii::VectorOperation::add);
      
    }
    

    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    unsigned int
    Model<dim, GroupSymm, ValueType>::solve ()
    {
      
      // Ultimately, a vector is required that stores not only the
      // elements of the solution for degrees of freedom the current
      // processor owns, but also all other locally relevant degrees
      // of freedom. The solver itself needs a vector that is uniquely
      // split between processors, without any overlap. Therefore
      // create a vector at the beginning of this function that has
      // these properties, use it to solve the linear system, and only
      // assign it to the vector we want at the very end. This last
      // step ensures that all ghost elements are also copied as
      // necessary.
      dealii::SolverControl solver_control (system_matrix.m (), 1e-09*system_rhs.l2_norm ());
#ifdef USE_PETSC
      dealii::PETScWrappers::MPI::Vector distributed_solution (locally_owned_dofs, mpi_communicator);
      dealii::PETScWrappers::PreconditionBlockJacobi preconditioner (system_matrix);
      dealii::PETScWrappers::SolverBicgstab solver (solver_control, mpi_communicator);
#else
      dealii::TrilinosWrappers::MPI::Vector distributed_solution (locally_owned_dofs, mpi_communicator);
      dealii::TrilinosWrappers::PreconditionBlockJacobi preconditioner;
      preconditioner.initialize (system_matrix);
      dealii::TrilinosWrappers::SolverBicgstab solver (solver_control);
#endif
      
      solver.solve (system_matrix, distributed_solution, system_rhs, 
		    preconditioner);
      
      constraints.distribute (distributed_solution);
      solution = distributed_solution;
      
      return solver_control.last_step ();
    }
    
    
    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    void
    Model<dim, GroupSymm, ValueType>::refine_grid ()
    {
      dealii::Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
      
      // Attempt an estimate based on Kelly's error estimate of each
      // component of the solution vector. @todo The Gauss quadrature
      // used here should be probably greater than that of the actual
      // solution vector.
      dealii::KellyErrorEstimator<dim>::estimate (dof_handler,
						  dealii::QGauss<dim-1> (9),
						  typename dealii::FunctionMap<dim>::type (),
						  solution, estimated_error_per_cell);
      
      // Setup grid refinement by fixed numbers.
      dealii::parallel::distributed::GridRefinement::
	refine_and_coarsen_fixed_number (triangulation,
					 estimated_error_per_cell,
					 0.125, 0.00);
      
      // Actually refine the grid.
      triangulation.execute_coarsening_and_refinement ();
    }
    
    
    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    void
    Model<dim, GroupSymm, ValueType>::output_material_id (const unsigned int cycle) const
    {
      
      // This is a speedy relatively operation, so we can get away
      // with running it on one processor only. In fact, we should
      // attach a different dof handler to this since the current one
      // has four components but only one is needed for material_id.
      dealii::DataOut<dim> data_out;
      data_out.attach_dof_handler (dof_handler);
      
#ifdef USE_PETSC
      dealii::PETScWrappers::Vector projected_material_id (dof_handler.n_dofs ());
#else
      dealii::TrilinosWrappers::Vector projected_material_id (dof_handler.n_dofs ());
#endif
      
      // nil::GeometryFunction::SquarePyramid<dim, ValueType> geometry (10., 5., 5.);
      nil::GeometryFunction::HalfHyperBall<dim, ValueType> geometry (5., dealii::Point<dim, ValueType> ());
      dealii::VectorTools::interpolate (dof_handler, geometry, 
					projected_material_id);
      
      const std::string filename 
	= ("projected_material_id-" +
	   dealii::Utilities::int_to_string (cycle, 4) +
	   "." +
	   dealii::Utilities::int_to_string (triangulation.locally_owned_subdomain(), 4) +
	   ".vtu");
      
      data_out.add_data_vector (projected_material_id, "projected_material_id");
      data_out.build_patches ();
      
      std::ofstream output (filename.c_str()); 
      data_out.write_vtu (output);
      
      // For calculations in parallel it is convenient to have a
      // master record written. This allows paravew (and friends) to
      // read the entire structure by binding blocks (*.vtu) into a
      // single structure (*.visit). This only needs to be done by one
      // processor.
      if (dealii::Utilities::MPI::this_mpi_process (mpi_communicator) == 0)
	{
	  std::vector<std::string> filenames;
	  
	  // The filenames to bind into the master record should match
	  // the filenames used above.
	  for (unsigned int i=0; i<dealii::Utilities::MPI::n_mpi_processes (mpi_communicator); ++i)
	    {
	      filenames.push_back ("projected_material_id-" +
				   dealii::Utilities::int_to_string (cycle, 4) +
				   "." +
				   dealii::Utilities::int_to_string (i, 4) +
				   ".vtu");
	    }
	  
	  const std::string
	    pvtu_master_filename = ("projected_material_id-" +
				    dealii::Utilities::int_to_string (cycle, 4) +
				    ".pvtu");
	  
	  std::ofstream pvtu_master (pvtu_master_filename.c_str ());
	  data_out.write_pvtu_record (pvtu_master, filenames);
	}
      
    }
    
    
    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    void
    Model<dim, GroupSymm, ValueType>::output_results (const unsigned int cycle) const
    {
      Model::Postprocessor postprocessor (dealii::Utilities::MPI::this_mpi_process (mpi_communicator));
      
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
      
      // For calculations in parallel it is convenient to have a
      // master record written. This allows paravew (and friends) to
      // read the entire structure by binding blocks (*.vtu) into a
      // single structure (*.visit). This only needs to be done by one
      // processor.
      if (dealii::Utilities::MPI::this_mpi_process (mpi_communicator) == 0)
	{
	  std::vector<std::string> filenames;
	  
	  // The filenames to bind into the master record should match
	  // the filenames used above.
	  for (unsigned int i=0; i<dealii::Utilities::MPI::n_mpi_processes (mpi_communicator); ++i)
	    {
	      filenames.push_back ("solution-" +
				   dealii::Utilities::int_to_string (cycle, 4) +
				   "." +
				   dealii::Utilities::int_to_string (i, 4) +
				   ".vtu");
	    }
	  
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
    
    
    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    void 
    Model<dim, GroupSymm, ValueType>::run ()
    {
      
      // First find the parameters need for this calculation
      get_parameters ();
      
      // Make coefficient tensors.
      setup_coefficient_tensors ();
      
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
      
      // Then create the coarse grid
      make_coarse_grid (3);
      
      // Here comes the adaptive cycles
      const unsigned int n_cycles = 3;
      
      for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
	{
	  
	  pcout << std::endl
		<< "Cycle:                              "
		<< cycle
		<< std::endl;
	  
	  pcout << "   Number of active cells:          "
		<< triangulation.n_global_active_cells ()
		<< " (on "
		<< triangulation.n_levels ()
		<< " levels)"
		<< std::endl;
	  
	  // Distribute degrees of freedom and reinitialize matrices and
	  // vectors.
	  setup_system ();
	  
	  pcout << "   Number of degrees of freedom:    "
		<< dof_handler.n_dofs ()
		<< std::endl;
	  
	  // Assemble the matrix and rhs vector.
	  assemble_system ();
	  
	  pcout << "   Number of non-zero elements:     "
		<< system_matrix.n_nonzero_elements ()
		<< std::endl
		<< "   System rhs l2-norm:              "
		<< system_rhs.l2_norm ()
		<< std::endl;
	  
	  // Solve the problem.
	  const unsigned int n_iterations = solve ();
	  
	  pcout << "   Solver converged in:             "
		<< n_iterations << " iterations"
		<< std::endl;
	  
	  // Write derived quantities to files.
	  output_results (cycle);
	  output_material_id (cycle);
	  
	  // Finally use Kelly's error estimate to refine the grid.
	  refine_grid ();
	}
      
    }
    
  } // namespace Piezoelectric
  
} // namespace nil


// -------------- Explicit Instantiations -------------------------------

template class 
nil::Piezoelectric::Model<3, nil::GroupSymmetry::ZincBlende, double>;

template class 
nil::Piezoelectric::Model<3, nil::GroupSymmetry::Wurtzite, double>;
