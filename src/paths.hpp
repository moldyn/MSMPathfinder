#ifndef PATHS_H
#define PATHS_H

#include <iostream>
#include <fstream>
#include "main.hpp"
#include "pathfinder.hpp"
#include "transition_matrix.hpp"
#include "progressbar.hpp"

#include <string>
#include <vector>
#include <cmath>

//#include <exception>
//#include <omp.h>
void generate_pathfinders(const Mode mode,
                          Pathfinder_map &pathfinders,
                          const std::vector<std::size_t> &A,
                          const std::vector<std::size_t> &B,
                          const std::vector<std::size_t> &C,
                          const Transition_matrix &T,
                          const std::size_t steps,
                          const std::size_t total_steps,
                          const std::size_t iterations,
                          const double cut_off,
                          const std::size_t threshold,
                          const bool isWeight,
                          const bool normalize);
int run_paths(Mode mode,
              std::vector<std::size_t> A,
              std::vector<std::size_t> B,
              std::vector<std::size_t> C,
              double cut_off,
              std::size_t steps,
              std::size_t iterations,
              std::size_t total_steps,
              Transition_matrix T,
              std::string output,
              std::size_t threshold,
              int argc, char* argv[]);

#endif
