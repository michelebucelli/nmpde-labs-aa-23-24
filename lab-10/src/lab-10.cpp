#include "Poisson2D.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  Poisson2D problem_0(0);
  Poisson2D problem_1(1);

  problem_0.setup();
  problem_1.setup();

  std::cout << "Setup completed" << std::endl;

  const double       tolerance_increment = 1e-4;
  const unsigned int n_max_iter          = 100;

  double       solution_increment_norm = tolerance_increment + 1;
  unsigned int n_iter                  = 0;

  while (n_iter < n_max_iter && solution_increment_norm > tolerance_increment)
    {
    }

  return 0;
}