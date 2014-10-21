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

#include "six_band_hole_model.h"

namespace nil
{

  namespace SixBandHole
  {
    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    Model<dim, GroupSymm, ValueType>::Model (dealii::parallel::distributed::Triangulation<dim> &coarse_grid,
					     MPI_Comm                                           mpi_communicator)
      :
      mpi_comm (mpi_communicator),
      
      pcout (std::cout, (dealii::Utilities::MPI::this_mpi_process (mpi_communicator) == 0)),
      
      triangulation (&coarse_grid),
      
      fe_q (dealii::FE_Q<dim> (2), 6), /* six holes */
      
      dof_handler (*triangulation),
      
      n_material_ids (2)
    {}
    
    template <int dim, enum nil::GroupSymmetry GroupSymm, typename ValueType>
    Model<dim, GroupSymm, ValueType>::~Model ()
    {
      dof_handler.clear ();
    }
    
  } // namespace SixBandHole

} // namespace nil

// -------------- Explicit Instantiations -------------------------------

template class 
nil::SixBandHole::Model<3, nil::GroupSymmetry::ZincBlende, double>;

template class 
nil::SixBandHole::Model<3, nil::GroupSymmetry::Wurtzite, double>;
