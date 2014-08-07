
// deal.II headers
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
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
  
  // First is a list of functions that belong to this class (they are
  // explained later on).
  void setup_problem ();
  void get_parameters ();
  
  // Following that we have a list of the tensors that will be used in
  // this calculation. They are, first- and second-order piezoelectric
  // tensors, a first-order strain tensor, and a tensor of first-order
  // displacement
  nil::PiezoelectricTensor<1, ValueType> first_order_piezoelectric_tensor;
  nil::PiezoelectricTensor<2, ValueType> second_order_piezoelectric_tensor;
  
  dealii::Tensor<1, dim> first_order_displacement;
  dealii::Tensor<1, dim> first_order_strain;
  
  // Additionally, lists of coefficients are needed for those tensors
  // that are tensors of empirical moduli.
  std::vector<ValueType> first_order_piezoelectric_coefficients;
  std::vector<ValueType> second_order_piezoelectric_coefficients;

  // Then we need an object to hold various run-time parameters that
  // are specified in an "prm file".
  dealii::ParameterHandler parameters;
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
  first_order_piezoelectric_tensor.reinit ();
  second_order_piezoelectric_tensor.reinit ();

  // and then the strain tensor.
  first_order_strain.reinit ();
}


// The next step is to obtain a complete list of the coefficients
// needed for this calculation. First comes a declaration of the
// entries expected to be find in the parameter file and then they are
// read into the object parameters.
template <int dim, typename ValueType>
void 
Step0<dim, ValueType>::get_parameters ()
{
  // First declare the parameters that are expected to be found.
  parameters.declare_entry ("First-order piezoelectric constants",
			    "-0.230",
			    dealii::Patterns::List (dealii::Patterns::Double (), 1),
			    "A list of the first-order piezoelectric constants. " 
			    "Default is zinc-blende GaAs. "
			    "See: PRL 96, 187602 (2006). ");
  
  parameters.declare_entry ("Second-order piezoelectric constants",
			    "-0.439, -3.765, -0.492",
			    dealii::Patterns::List (dealii::Patterns::Double (), 1),
			    "A list of the second-order piezoelectric constants. "
			    "Default is zinc-blende GaAs. "
			    "See: PRL 96, 187602 (2006). ");

  parameters.read_input (command_line.get_prm_file ());

  first_order_piezoelectric_coefficients = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("First-order piezoelectric constants")));

  second_order_piezoelectric_coefficients = 
    dealii::Utilities::string_to_double
    (dealii::Utilities::split_string_list (parameters.get ("Second-order piezoelectric constants")));
}


template <int dim, typename ValueType>
void 
Step0<dim, ValueType>::run ()
{
  // First find the parameters need for this calculation
  get_parameters ();

  // First up, fill the piezoelectric tensors with coefficent
  // values. Starting with first-order coefficients...  
  first_order_piezoelectric_tensor
    .distribute_first_order_piezoelectric_coefficients (first_order_piezoelectric_coefficients);
  
  // and then second-order coefficients.
  second_order_piezoelectric_tensor
    .distribute_second_order_piezoelectric_coefficients (second_order_piezoelectric_coefficients);
  
  // Having done that now we want to start applying an incremental
  // strain pattern. This is done 
  
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
