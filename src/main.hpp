#include <iostream>
#include <fstream>
#include "pathfinder.hpp"
#include "transition_matrix.hpp"
#include "progressbar.hpp"

#include <string>
#include <vector>
#include <cmath>

//#include <exception>

#include <boost/program_options.hpp>
//#include <omp.h>

#ifndef MODE
#define MODE
enum Mode {PATHS, MCMC, MCMCSINGLE};
typedef std::map<std::size_t, double> Weights;


//! check if states are legit values and unique
bool check_states(const std::vector<std::size_t> A, const std::string str_A,
                  const std::vector<std::size_t> B, const std::string str_B,
                  const std::vector<std::size_t> C, const std::string str_C,
                  const std::size_t no_of_states);

//! check if states in A are legit
bool check_states_A(const std::vector<std::size_t> A, const std::string str_A,
                    const std::vector<std::size_t> B, const std::string str_B,
                    const std::vector<std::size_t> C, const std::string str_C,
                    const std::size_t no_of_states);


#endif
