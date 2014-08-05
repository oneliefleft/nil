
// deal.II headers
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>

// Library-based headers.
#include "piezoelectric_tensor.h"

// Application-base headers.
#include "command_line.h"

#include <fstream>
#include <iostream>


template <int dim, typename ValueType = double>
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
  
  // First is a list of functions that belong to this class.
  void setup_problem ();
  
  // Following that we have a list of the tensors that will be used in
  // this calculation. They are, first- and second-order piezoelectric
  // tensors, a first-order strain tensor, and a tensor of first-order
  // displacement
  nil::PiezoelectricTensor<3, 1, ValueType> first_order_piezoelectric;
  nil::PiezoelectricTensor<3, 2, ValueType> second_order_piezoelectric;
  
  dealii::Tensor<1, dim> first_order_displacement;
  dealii::Tensor<1, dim> first_order_strain;
  
  // Additionally, lists of coefficients are needed for those tensors
  // that are tensors of empirical moduli.
  std::list<ValueType> first_order_piezoelectric_coefficients;
  std::list<ValueType> second_order_piezoelectric_coefficients;

};

// The constructor is typically borning...
template <int dim, typename ValueType>
Step0<dim, ValueType>::Step0 ()
{}

// as is the destructor.
template <int dim, typename ValueType>
Step0<dim, ValueType>::~Step0 ()
{}


template <int dim, typename ValueType>
void 
Step0<dim, ValueType>::setup_problem ()
{
  // Initialise the first- and second-order piezoelectric tensors...
  first_order_piezoelectric.reinit ();
  second_order_piezoelectric.reinit ();

  // and then the strain tensor.
  first_order_strain.reinit ();
}

template <int dim, typename ValueType>
void 
Step0<dim, ValueType>::run ()
{
  // First up, fill the piezoelectric tensors with coefficent
  // values. Starting with first-order coefficients...  
  first_order_piezoelectric
    .distribute_first_order_piezoelectric_coefficients (first_order_piezoelectric_coefficients);
  
  // and then second-order coefficients.
  second_order_piezoelectric
    .distribute_second_order_piezoelectric_coefficients (second_order_piezoelectric_coefficients);
  
  // Having done that now we want to start applying an incremental
  // strain pattern
  
}


int main (int argc, char **argv)
{

  try
    {
      // Set deal.II's logger depth.
      dealii::deallog.depth_console (0);

      // Initialise the problem.
      Step0<3, double> problem;

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
