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
                          const std::vector<State> &A,
                          const std::vector<State> &B,
                          const std::vector<State> &C,
                          const Transition_matrix &T,
                          const std::size_t steps,
                          const std::size_t total_steps,
                          const std::size_t iterations,
                          const double cut_off,
                          const std::size_t threshold,
                          const bool isWeight,
                          const bool normalize);
int run_paths(Mode mode,
              std::vector<State> A,
              std::vector<State> B,
              std::vector<State> C,
              double cut_off,
              std::size_t steps,
              std::size_t iterations,
              std::size_t total_steps,
              Transition_matrix T,
              std::string output,
              std::size_t threshold,
              int argc, char* argv[]);

#endif
